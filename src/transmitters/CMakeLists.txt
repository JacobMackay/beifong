set(MTS_PLUGIN_PREFIX "transmitters")

add_plugin(wignertransmitter   wignertransmitter.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
