function result = run_performance_tests()
    import matlab.unittest.TestSuite

    performace_tests = TestSuite.fromClass(?PerformanceTests);
    
    fullSuite = performace_tests;
    result = run(fullSuite);
end