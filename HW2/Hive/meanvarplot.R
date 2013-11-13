#STA 250 HOMEWORK 2
#PROBLEM 3
#PLOT WITHIN GROUP MEANS AGAINST WITHIN GROUP VARIANCES

dat = read.table("~/Desktop/SkyDrive/STA 250/Stuff/HW2/Hive/group_results.txt",
                 header = FALSE,sep = ",")

png("~/Desktop/SkyDrive/STA 250/Stuff/HW2/Hive/mean_var.png",
    width=800, height = 600)
plot(dat[,2]~dat[,1], type = "p",
     main = "Within-group variances vs within-group means",
     xlab ="means", ylab = "variances")
dev.off()




