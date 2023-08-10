function result = run_all_tests()
    import matlab.unittest.TestSuite

    glodap_tests = TestSuite.fromClass(?GlodapTests);
    grid_tests = TestSuite.fromClass(?GridTests);
    random_tests = TestSuite.fromClass(?RandomTests);
    performace_tests = TestSuite.fromClass(?PerformanceTests);
    chain_test = TestSuite.fromClass(?ChainTests);
    
    fullSuite = [glodap_tests,grid_tests,random_tests,performace_tests,chain_test];
    result = run(fullSuite);
end