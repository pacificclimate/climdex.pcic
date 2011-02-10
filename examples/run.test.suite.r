library(RUnit)

test.suite <- defineTestSuite("pcic_climdex",
                                 dirs = c("./"),
                                 testFileRegexp = "^test.+r$",
                                 testFuncRegexp = "^test.+")

test.result   <- runTestSuite(test.suite, useOwnErrorHandler=F)
printTextProtocol(test.result)
