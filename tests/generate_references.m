% Glodap
glodap_tests = GlodapTests();

glodap_tests.generate_glodap_results
glodap_tests.generate_glodap_data_subset();
glodap_tests.generate_glodap_subset_results();

% Grid
grid_tests = GridTests();

grid_tests.generate_temperature_salinity_pressure_reference();
grid_tests.generate_dic_alkalinity_reference();

% Random
random_tests = RandomTests();

random_tests.generate_temperature_salinity_pressure_reference();
random_tests.generate_dic_alkalinity_reference();