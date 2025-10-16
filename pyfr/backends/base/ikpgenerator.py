import re

from pyfr.backends.base.makoutil import C_TYPE_QUALIFIERS


class HoistVar:
    """Represents a local variable that needs hoisting to block scope."""
    def __init__(self, dtype, name, dimstr, qual=''):
        self.dtype = dtype
        self.name = name
        self.cdimstr = dimstr
        self.qual = qual  # e.g., 'const', '__constant__ const', etc.

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
    # Match with optional qualifiers: [qualifiers] type name[dims];
    # Captures: (qualifiers, type, name, dims)
    _LOCAL_ARRAY_PATTERN = (
        rf'((?:(?:{"|".join(C_TYPE_QUALIFIERS)})\s+)*)'  # Qualifiers (capture group 1)
        r'([A-Za-z_]\w*)\s+'                              # Type (capture group 2)
        r'(\w+)'                                          # Name (capture group 3)
        r'((?:\[\d+\])+)'                                # Dimensions (capture group 4)
        r'\s*(?:=\s*\{[^}]*\})?\s*;'                     # Optional initializer
    )
    # Match with optional qualifiers: [qualifiers] type name;
    # Captures: (qualifiers, type, name)
    _LOCAL_SCALAR_PATTERN = (
        rf'((?:(?:{"|".join(C_TYPE_QUALIFIERS)})\s+)*)'  # Qualifiers (capture group 1)
        r'([A-Za-z_]\w*)\s+'                              # Type (capture group 2)
        r'(\w+)'                                          # Name (capture group 3)
        r'(?!\s*\[)'                                      # Not followed by [
        r'\s*(?:=\s*[^;]+)?;'                            # Optional initializer
    )

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

    def _get_kernel_level_code(self, content):
        """
        Extract only section-level code (outside all {...} blocks).

        Returns code that appears before the first { and after the last },
        which represents declarations at section scope.
        """
        # Find first { and last }
        first = content.find('{')
        last = content.rfind('}')

        if first == -1:
            # No braces, entire content is section-level
            return content

        # Section-level code is before first { and after last }
        before = content[:first] if first != -1 else ''
        after = content[last + 1:] if last != -1 else ''

        return before + '\n' + after

    def _find_hvars(self, secs):
        """
        Find local variables that need hoisting to block scope.

        Variables are hoisted if:
        - Declared at section-level (outside all {...} blocks) in a looped section
        - Used in any downstream section

        Raises:
            ValueError: If variable declaration intent is ambiguous

        Returns:
            list: HoistVar objects for variables needing hoisting
        """
        hvars = []
        seen = {}  # Track where each name was declared

        # Process each looped section
        for idx, (stype, content) in enumerate(secs):
            if stype == 'noloop':
                continue

            # Get only kernel-level code (outside scoped nests)
            klevel = self._get_kernel_level_code(content)

            # Collect all kernel-level declarations (arrays and scalars)
            decls = []
            for match in re.finditer(self._LOCAL_ARRAY_PATTERN, klevel):
                quals, dtype, name, dimstr = match.groups()
                decls.append((quals, dtype, name, dimstr))

            for match in re.finditer(self._LOCAL_SCALAR_PATTERN, klevel):
                quals, dtype, name = match.groups()
                decls.append((quals, dtype, name, ''))

            # Process all declarations
            for quals, dtype, name, dimstr in decls:
                # Check for duplicate declarations
                if name in seen:
                    raise ValueError(
                        f'IKP: Variable "{name}" declared in multiple sections'
                        f' (sections {seen[name]} and {idx}). '
                        f'Hoisting intent is ambiguous.'
                    )
                seen[name] = idx

                # Check if used in any downstream section (full content)
                for dsidx in range(idx + 1, len(secs)):
                    if re.search(rf'\b{name}\b', secs[dsidx][1]):
                        hvars.append(HoistVar(dtype, name, dimstr, quals))
                        break

        return hvars
