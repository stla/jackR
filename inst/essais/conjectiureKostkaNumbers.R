library(jack)
library(gmp)

n <- 5 # weight of the partitions

partitions_n     <- partitions::parts(n)
dualPartitions_n <- partitions::conjugate(partitions_n)

partitions_n_asStrings     <- apply(partitions_n, 2L, toString)
dualPartitions_n_asStrings <- apply(dualPartitions_n, 2L, toString)

KJN_three    <- KostkaJackNumbers(n, alpha = "3")
KJN_oneThird <- KostkaJackNumbers(n, alpha = "1/3")

rownames(KJN_oneThird) <- colnames(KJN_oneThird) <- partitions_n_asStrings
KJN_oneThird <- KJN_oneThird[dualPartitions_n_asStrings, dualPartitions_n_asStrings]

t(as.bigq(t(KJN_three) %*% as.bigq(KJN_oneThird))

KostkaJackNumbers(n, alpha = "0")
