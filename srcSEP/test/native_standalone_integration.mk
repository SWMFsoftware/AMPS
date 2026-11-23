# Use the ordinary configured CPU field-line build's own compiler/link flags.
# No reconfiguration, physical timestep override or MPI dependency substitute.
include Makefile
.PHONY: native-standalone-integration
native-standalone-integration:
	$(CC) $(SEARCH_C) -std=c++17 -O0 "$${AMPS_INTEGRATION_REPO}/srcSEP/test/native_standalone_integration.cpp" $(WSD)/libAMPS.a $(CPPLIB) $(AMPSLINKLIB) $(EXTRALINKEROPTIONS) -o "$${AMPS_INTEGRATION_OUTPUT}"

# The same archive, PIC ABI and real MPI library as the production consumer.
# Only operator fixtures live in this executable; no alternative physics or
# replacement MPI implementation is linked to make a qualification gate pass.
.PHONY: native-checkpoint-operators
native-checkpoint-operators:
	$(CC) $(SEARCH_C) -std=c++17 -O0 "$${AMPS_INTEGRATION_REPO}/srcSEP/test/native_checkpoint_operators.cpp" $(WSD)/libAMPS.a $(CPPLIB) $(AMPSLINKLIB) $(EXTRALINKEROPTIONS) -o "$${AMPS_INTEGRATION_OUTPUT}"
