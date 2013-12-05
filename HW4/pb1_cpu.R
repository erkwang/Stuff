#STA 250 HOMEWORK 4
#YICHUAN WANG

#PROBLEM 1

#TRUNCATED NORM WITH CPU

library(msm)

cpu_time = list()
cpu_values = list()
for(k in 1:8){
N = 10^k
cpu_time[[k]] = system.time({
  cpu_values[[k]] = rtnorm(N, mean=2, sd = 1, lower = 0, upper=1.5)
})
}
cpu_results = list(cpu_values, cpu_time)
save(cpu_results, file = "~/Desktop/cpu_results.RDA")