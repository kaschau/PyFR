import re


class HoistVar:
    """Represents a local variable that needs hoisting to block scope."""
    def __init__(self, dtype, name, dimstr):
        self.dtype = dtype
        self.name = name
        self.cdimstr = dimstr

        # Parse dimensions from string like "[3][4]"
        dimsptn = r'(?<=\[)\d+(?=\])'
        self.cdims = [int(d) for d in re.findall(dimsptn, dimstr)]
        self.ncdim = len(self.cdims)

        # Classify variable type
        self.isscalar = self.ncdim == 0
        self.isarray = self.ncdim > 0


class IKPKernelGeneratorMixin:
    """
    Mixin class providing common IKP transformation utilities.

    This class handles:
    - Pattern matching for constant and local array declarations
    - Tracking which arrays need IKP transformation
    - Common regex patterns

    Backend-specific subclasses implement:
    - How to transform declarations (stack vs shared memory)
    - How to transform references (ELEM_IDX, threadIdx.x, etc.)
    - How to wrap body (sequential loop vs thread parallelism)
    """

    # Regex patterns for declaration matching (common to all backends)
    # Match any type: int name[5]; fpdtype_t arr[3][4]; double mat[2][3][4]; etc.
    _LOCAL_ARRAY_PATTERN = r'\s*([A-Za-z_]\w*)\s+(\w+)((?:\[\d+\])+)\s*(?:=\s*\{[^}]*\})?\s*;'
    # Match any type: int i; fpdtype_t x = expr; size_t n;
    _LOCAL_SCALAR_PATTERN = r'\s*([A-Za-z_]\w*)\s+(\w+)(?!\s*\[)\s*(?:=\s*[^;]+)?;'

    def _split_ikp_sections(self, body):
        """
        Split IKP body into sections based on loop requirements.

        Sections are delimited by **IKP_SECTION_START/END markers.
        Each section is classified as 'looped' or 'noloop' based on
        presence of // NO LOOP marker.

        Returns:
            list of (section_type, content) tuples where section_type is:
            - 'looped': needs per-element loop wrapping
            - 'noloop': batched operation, no loop needed
        """
        # Extract all sections between START/END markers
        pattern = r'\*\*IKP_SECTION_START\s*(.*?)\s*\*\*IKP_SECTION_END'
        matches = re.finditer(pattern, body, re.S)

        sections = []
        for match in matches:
            content = match[1].strip()

            # Skip empty sections
            if not content:
                continue

            # Check for NO LOOP marker
            # HACK
            stype = 'noloop' if '// NO LOOP' in content else 'looped'
            sections.append((stype, content))

        return sections

    def _find_hvars(self, secs):
        """
        Find local variables that need hoisting to block scope.

        Variables are hoisted if declared in one section and used in another.

        Returns:
            list: HoistVar objects for variables needing hoisting
        """
        hvars = []

        # Process each looped section
        for idx, (stype, content) in enumerate(secs):
            if stype == 'noloop':
                continue

            # Find arrays declared in this section
            for match in re.finditer(self._LOCAL_ARRAY_PATTERN, content):
                dtype, name, dimstr = match.groups()

                # Check if used in any downstream section
                for dsidx in range(idx + 1, len(secs)):
                    if re.search(rf'\b{name}\b', secs[dsidx][1]):
                        hvars.append(HoistVar(dtype, name, dimstr))
                        break

            # Find scalars declared in this section
            for match in re.finditer(self._LOCAL_SCALAR_PATTERN, content):
                dtype, name = match.groups()

                # Check if used in any downstream section
                for dsidx in range(idx + 1, len(secs)):
                    if re.search(rf'\b{name}\b', secs[dsidx][1]):
                        hvars.append(HoistVar(dtype, name, ''))
                        break

        return hvars
