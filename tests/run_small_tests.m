function result = run_small_tests()
    import matlab.unittest.TestSuite

    glodap_test_1 = TestSuite.fromMethod(?GlodapTests,"compare_glodap_subset");
    grid_tests = TestSuite.fromClass(?GridTests);
    random_tests = TestSuite.fromClass(?RandomTests);
    chain_test = TestSuite.fromClass(?ChainTests);
    equilibrium_constant_tests = TestSuite.fromClass(?EquilibriumConstantsTests);
    length_tests = TestSuite.fromClass(?LengthTests);
    
    fullSuite = [glodap_test_1,...
                 grid_tests,...
                 random_tests,...
                 chain_test,...
                 equilibrium_constant_tests,...
                 length_tests];
    result = run(fullSuite);
end