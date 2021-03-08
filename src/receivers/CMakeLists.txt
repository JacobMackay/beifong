set(MTS_PLUGIN_PREFIX "receivers")

# add_plugin(perspective     perspective.cpp)
# add_plugin(radiancemeter   radiancemeter.cpp)
# add_plugin(thinlens        thinlens.cpp)
# add_plugin(irradiancemeter irradiancemeter.cpp)
# add_plugin(fluxmeter fluxmeter.cpp)

add_plugin(omnidirectional   omnidirectional.cpp)
add_plugin(wignerreceiver   wignerreceiver.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
