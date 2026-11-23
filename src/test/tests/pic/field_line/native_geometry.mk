# Include the configured build's ordinary compiler, preprocessor and linker
# settings. This fragment never reconfigures the checkout or changes its dt.
include Makefile

.PHONY: native-field-line-fixture native-field-line-fixture-serial
native-field-line-fixture:
	$(CC) $(SEARCH_C) -std=c++17 -O0 "$${AMPS_NATIVE_TEST_REPO}/src/test/tests/pic/field_line/native_geometry.cpp" $(WSD)/libAMPS.a $(CPPLIB) $(AMPSLINKLIB) $(EXTRALINKEROPTIONS) -o "$${AMPS_NATIVE_TEST_OUTPUT}"

# The adapter executable is deliberately a different target and test identity.
# A successful serial dependency test never supplies an MPI qualification gate.
native-field-line-fixture-serial:
	$(CC) $(SEARCH_C) -std=c++17 -O0 "$${AMPS_NATIVE_TEST_REPO}/src/test/tests/pic/field_line/native_geometry.cpp" "$${AMPS_NATIVE_TEST_REPO}/src/test/tests/pic/field_line/serial_mpi_dependency.cpp" $(WSD)/libAMPS.a $(CPPLIB) $(AMPSLINKLIB) $(EXTRALINKEROPTIONS) -o "$${AMPS_NATIVE_TEST_OUTPUT}"
