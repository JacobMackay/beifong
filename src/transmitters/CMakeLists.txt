set(MTS_PLUGIN_PREFIX "transmitters")

add_plugin(wignertransmitter   wignertransmitter.cpp)
add_plugin(phasedtransmitter   phasedtransmitter.cpp)
add_plugin(areatransmitter   areatransmitter.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
