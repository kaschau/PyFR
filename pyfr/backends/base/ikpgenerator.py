import re


class IKPLocalVar:
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
    _LOCAL_SCALAR_PATTERN = r'\s*([A-Za-z_]\w*)\s+(\w+)\s*(?:=\s*[^;]+)?;'

    def _ikp_split_into_sections(self, body):
        """
        Split IKP body into alternating prep/interruption sections.

        Returns:
            list of (section_type, content) tuples
            section_type is 'prep' or 'interruption'
        """
        sections = []
        remaining = body

        while True:
            # Find next interruption
            match = re.search(
                r'(.*?)// PYFR_IKP_INTERRUPTION_START\n(.*?)// PYFR_IKP_INTERRUPTION_END',
                remaining, flags=re.DOTALL
            )

            if not match:
                # No more interruptions, rest is final prep/proc section
                if remaining.strip():
                    sections.append(('prep', remaining))
                break

            # Extract prep section before interruption
            prep_section = match.group(1)
            if prep_section.strip():
                sections.append(('prep', prep_section))

            # Extract interruption section
            interruption_section = match.group(2)
            sections.append(('interruption', interruption_section))

            # Continue with remainder
            remaining = remaining[match.end():]

        return sections

    def _ikp_analyze_variable_usage(self, sections):
        """
        Analyze which variables are used in which sections.

        Returns:
            dict: {var_name: {'localvar': IKPLocalVar instance,
                              'declared_in': section_index,
                              'used_in': set of section_indices}}
        """
        variables = {}

        for idx, (section_type, content) in enumerate(sections):
            # Skip interruptions - they don't have local declarations
            if section_type == 'interruption':
                continue

            # Find local array declarations (1D or 2D)
            for match in re.finditer(self._LOCAL_ARRAY_PATTERN, content):
                dtype, name, dimstr = match.groups()
                if name not in variables:
                    variables[name] = {
                        'localvar': IKPLocalVar(dtype, name, dimstr),
                        'declared_in': idx,
                        'used_in': set()
                    }

            # Find local scalar declarations (but filter out arrays)
            content_no_arrays = re.sub(self._LOCAL_ARRAY_PATTERN, '', content)

            for match in re.finditer(self._LOCAL_SCALAR_PATTERN, content_no_arrays):
                dtype, name = match.groups()
                if name not in variables:
                    variables[name] = {
                        'localvar': IKPLocalVar(dtype, name, ''),  # Empty dimstr for scalars
                        'declared_in': idx,
                        'used_in': set()
                    }

        # Now find usages of each variable across all sections
        for var_name in variables.keys():
            for idx, (section_type, content) in enumerate(sections):
                # Check if variable is used in this section (look for var_name as a word)
                if re.search(r'\b' + re.escape(var_name) + r'\b', content):
                    variables[var_name]['used_in'].add(idx)

        return variables

    # Backend-specific methods (must be implemented by subclasses)

    def _ikp_transform_local_decl(self, localvar):
        """
        Transform local array declaration for IKP.

        OpenMP: arr[N] -> arr[BLK_SZ*N] (flat stack allocation)
        GPU: arr[N] -> __shared__ arr[BLOCK_DIM*N] (flat shared memory)
        """
        raise NotImplementedError("Backend must implement _ikp_transform_local_decl")

    def _ikp_transform_array_ref(self, localvar, body):
        """
        Transform array references in IKP sections.

        OpenMP: arr[i] -> arr[i*BLK_SZ + X_IDX] (flat indexing)
        GPU: arr[i] -> arr[i*BLOCK_DIM + threadIdx.x] (flat indexing)
        """
        raise NotImplementedError("Backend must implement _ikp_transform_array_ref")

    def _ikp_transform_kernel_args(self, body):
        """
        Transform kernel argument references in IKP sections.

        OpenMP: X_IDX_AOSOA(...) -> ELEM_IDX
        GPU: Similar but might use threadIdx.x directly
        """
        raise NotImplementedError("Backend must implement _ikp_transform_kernel_args")

    def _ikp_wrap_body(self, body, nelem_expr):
        """
        Wrap transformed body in backend-specific parallelization.

        OpenMP: for (_elem = 0; _elem < nelem_expr; _elem++) { body }
        GPU: Just body (each thread is an element)
        """
        raise NotImplementedError("Backend must implement _ikp_wrap_body")
