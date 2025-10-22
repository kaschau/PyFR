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

    def _remove_hoisted_decls(self, body, hvars):
        """
        Remove hoisted variable declarations from body.

        Handles comma-separated declarations, qualifiers, arrays, scalars,
        and optional initializers. For scalars with initializers, preserves
        the assignment.

        Returns:
            str: Body with declarations removed
        """
        for hvar in hvars:
            if hvar.isarray:
                # For arrays, we must match the complete declaration with dimensions
                # This handles: dtype name[dims] [= {...}] [, other_vars...];
                # We remove just this variable from the declaration list

                # Match the array variable with optional initializer
                var_pattern = rf'{hvar.name}(?:\[\d+\])+\s*(?:=\s*\{{[^}}]*\}})?'

                # Try to match as part of comma-separated list
                # Case 1: var, rest...  (remove var and comma)
                pattern = rf'{hvar.qual}\s+{hvar.dtype}\s+{var_pattern}\s*,\s*'
                body = re.sub(pattern, f'{hvar.qual} {hvar.dtype} ', body)

                # Case 2: ..., var, rest...  (remove comma and var)
                pattern = rf',\s*{var_pattern}\s*(?=,|;)'
                body = re.sub(pattern, '', body)

                # Case 3: dtype var;  (only variable in declaration)
                pattern = rf'\s*{hvar.qual}\s+{hvar.dtype}\s+{var_pattern}\s*;'
                body = re.sub(pattern, '', body)
            else:
                # For scalars, preserve initializers if they exist
                # Case 1: dtype x = expr, rest...  (convert to x = expr;, keep rest as new declaration)
                pattern = rf'({hvar.qual}\s+{hvar.dtype}\s+){hvar.name}\s*=\s*([^,;]+)\s*,\s*'
                body = re.sub(pattern, rf'{hvar.name} = \2;\n\1', body)

                # Case 2: dtype x = expr;  (convert to x = expr;)
                pattern = rf'(\s*){hvar.qual}\s+{hvar.dtype}\s+{hvar.name}\s*=\s*([^;]+);'
                body = re.sub(pattern, rf'\1{hvar.name} = \2;', body)

                # Case 3: ..., x = expr, ...  (keep x = expr;, handle rest)
                pattern = rf',\s*{hvar.name}\s*=\s*([^,;]+)\s*(?=,|;)'
                body = re.sub(pattern, rf';\n{hvar.name} = \1', body)

                # Case 4: dtype x, rest...  (remove x, keep rest)
                pattern = rf'{hvar.qual}\s+{hvar.dtype}\s+{hvar.name}\s*,\s*'
                body = re.sub(pattern, f'{hvar.qual} {hvar.dtype} ', body)

                # Case 5: ..., x, rest...  (remove comma and x)
                pattern = rf',\s*{hvar.name}\s*(?=,|;)'
                body = re.sub(pattern, '', body)

                # Case 6: dtype x;  (remove entire declaration)
                pattern = rf'\s*{hvar.qual}\s+{hvar.dtype}\s+{hvar.name}\s*;'
                body = re.sub(pattern, '', body)

        return body


class GPUIKPKernelGeneratorMixin(IKPKernelGeneratorMixin):
    """
    GPU-specific IKP transformation logic (CUDA/HIP).

    Thread cooperation model (threadIdx.y cooperation):
    - Each X-thread (threadIdx.x) processes one element
    - Y-threads (threadIdx.y) cooperate within the element
    - Variables marked in **SHARED[...] are converted to shared memory
    - Uses block2d configuration: (32, 8, 1) → 32 elements, 8 cooperative threads
    """

    def _extract_shared_vars(self, body):
        """
        Extract variable names from **SHARED[...] markers.

        Returns:
            set: Variable names that need shared memory
        """
        names = set()
        for match in re.finditer(r'\*\*SHARED\[([^\]]+)\]', body):
            vstr = match[1]
            vnames = [v.strip() for v in vstr.split(',') if v.strip()]
            names.update(vnames)
        return names

    def _find_shared_var_decls(self, body, names):
        """
        Find declarations of variables marked for shared memory.

        Returns:
            list: HoistVar objects for shared variables
        """
        svars = []

        # Find array declarations
        for match in re.finditer(self._LOCAL_ARRAY_PATTERN, body):
            quals, dtype, name, dimstr = match.groups()
            if name in names:
                svars.append(HoistVar(dtype, name, dimstr, quals))

        # Find scalar declarations
        for match in re.finditer(self._LOCAL_SCALAR_PATTERN, body):
            quals, dtype, name = match.groups()
            if name in names:
                svars.append(HoistVar(dtype, name, '', quals))

        return svars

    def _decl_shared(self, svars):
        """
        Generate __shared__ declarations for shared variables.

        With threadIdx.y cooperation model:
        - Scalar: x -> __shared__ dtype x[blockDim.x];
        - 1D array: arr[N] -> __shared__ dtype arr[N*blockDim.x];
        - 2D array: arr[M][N] -> __shared__ dtype arr[M*blockDim.x][N];
        - ND array: arr[...] -> __shared__ dtype arr[first_dim*blockDim.x][remaining dims];

        Each X-thread (element) gets its own block of memory.

        Returns:
            list: Shared memory declaration strings
        """
        decls = []
        nelem_per_block = self.block2d[0]

        for svar in svars:
            if svar.isscalar:
                # Scalar: __shared__ dtype name[NELEM];
                decls.append(f'{self._shared_prfx} {svar.dtype} {svar.name}[{nelem_per_block}];')
            else:
                # Array: flatten first dimension with group dimension
                ldim = svar.cdims[0]
                tdims = ''.join(f'[{d}]' for d in svar.cdims[1:])
                decls.append(f'{self._shared_prfx} {svar.dtype} {svar.name}[{ldim}*{nelem_per_block}]{tdims};')

        return decls

    def _transform_shared_refs(self, body, svars):
        """
        Add thread index to all references.

        With threadIdx.y cooperation model:
        - Scalars: x -> x[_lid[0]]
        - 1D: arr[i] -> arr[i + _lid[0] * array_size]
        - 2D: arr[i][j] -> arr[i + _lid[0] * first_dim][j]
        - 3D: arr[i][j][k] -> arr[i + _lid[0] * first_dim][j][k]

        Each X-thread owns a contiguous block:
        - _lid[0]=0: arr[0..N-1]
        - _lid[0]=1: arr[N..2N-1]
        - etc.

        Returns:
            str: Body with transformed references
        """
        eidx = self._lid[0]  # Use backend-specific thread index

        for svar in svars:
            if svar.isscalar:
                # Scalar: name -> name[threadIdx.x]
                pattern = rf'\b{svar.name}\b'
                replacement = f'{svar.name}[{eidx}]'
                body = re.sub(pattern, replacement, body)
            else:
                # Array: build pattern matching all dimensions
                pattern = rf'\b{svar.name}' + r'\[([^\]]+)\]' * svar.ncdim

                # Build replacement: arr[i + threadIdx.x * array_size][j][k]...
                size = svar.cdims[0]
                replacement = rf'{svar.name}[(\1) + {eidx} * {size}]'
                for i in range(2, svar.ncdim + 1):
                    replacement += rf'[\{i}]'

                body = re.sub(pattern, replacement, body)

        return body

    def _ikp_render_body_preamble(self, body, preamble):
        """
        Apply GPU IKP transformations.

        Pipeline:
        1. Extract shared variable names from **SHARED[...] markers
        2. Strip all IKP markers (merge sections)
        3. Find declarations of shared variables
        4. Generate __shared__ declarations
        5. Remove original declarations from body
        6. Transform all references to add group index

        Returns:
            tuple: (transformed_body, preamble)
        """
        # Extract variables marked for shared memory (before stripping markers)
        snames = self._extract_shared_vars(body)

        # Strip all IKP markers to merge sections into one body
        body = re.sub(r'\*\*IKP_SECTION_START', '', body)
        body = re.sub(r'\*\*IKP_SECTION_END', '', body)
        body = re.sub(r'\*\*SHARED\[[^\]]+\]', '', body)

        # Find declarations of shared variables (in clean body)
        svars = self._find_shared_var_decls(body, snames)

        # Generate __shared__ declarations
        sdecls = self._decl_shared(svars)

        # Remove original declarations from body
        body = self._remove_hoisted_decls(body, svars)

        # Transform references to add group index
        body = self._transform_shared_refs(body, svars)

        # Add declarations to preamble
        if sdecls:
            preamble = '\n'.join(sdecls) + '\n' + preamble

        return body, preamble
