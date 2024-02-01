function X = doDataDenormalization(Xnorm, mean,std)

X = Xnorm*std+mean;