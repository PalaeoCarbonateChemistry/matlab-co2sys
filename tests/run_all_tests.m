function result = run_all_tests()
    import matlab.unittest.TestSuite

    glodap_tests = TestSuite.fromClass(?GlodapTests);
    grid_tests = TestSuite.fromClass(?GridTests);
    random_tests = TestSuite.fromClass(?RandomTests);
    
    fullSuite = [glodap_tests,grid_tests,random_tests];
    result = run(fullSuite);
end