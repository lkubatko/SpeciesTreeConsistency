source("comp_lik.R")

# Starting values for the branch lengths and theta
y = c(1.0,2.0,3.0,0.001) 

# The site pattern counts used as input as ordered as follows:
# xxxx
# xxxy
# xxyx
# xxyy
# zzxy
# xyxx
# xyxy
# zxzy
# yxxy
# yxxx
# xzzy
# zxyz
# xzyz
# xyzz
# xyzw

# There are 15 rooted trees on 4 tips, so we must consider 15 species trees

# The first input file, temp1.ssa, is ordered as
# Species1
# Species2
# Species3
# Species4

# Tree 1: ((1,2),(3,4))

data1 = matrix(scan("temp1.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = (data1[,2]+data1[,3])
n3 = (data1[,6]+data1[,10])
n4 = (data1[,7]+data1[,9])
n5 = data1[,4]
n6 = data1[,5]
n7 = data1[,14]
n8 = (data1[,8]+data1[,11]+data1[,12]+data1[,13])
n9 = data1[,15]

optim1 = optim(y,GetLikSymm,control=list(maxit=1500))
lik1 = optim1$value
conv1 = optim1$convergence
par1 = optim1$par

# Tree 2: (1,(2,(3,4)))

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]

optim2 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik2 = optim2$value
conv2 = optim2$convergence
par2 = optim2$par

# The second input file, temp2.ssa, is ordered as
# Species1
# Species3
# Species2
# Species4 

# Tree3: ((1,3),(2,4))

data1 = matrix(scan("temp2.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = (data1[,2]+data1[,3])
n3 = (data1[,6]+data1[,10])
n4 = (data1[,7]+data1[,9])
n5 = data1[,4]
n6 = data1[,5]
n7 = data1[,14]
n8 = (data1[,8]+data1[,11]+data1[,12]+data1[,13])
n9 = data1[,15]   

optim3 = optim(y,GetLikSymm,control=list(maxit=1500))
lik3 = optim3$value
conv3 = optim3$convergence
par3 = optim3$par

# Tree 4: (1,(3,(2,4)))

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim4 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik4 = optim4$value
conv4 = optim4$convergence
par4 = optim4$par

# The third input file, temp3.ssa, is ordered as
# Species1
# Species4
# Species2
# Species3     

# Tree 5: ((1,4),(2,3))

data1 = matrix(scan("temp3.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = (data1[,2]+data1[,3])
n3 = (data1[,6]+data1[,10])
n4 = (data1[,7]+data1[,9])
n5 = data1[,4]
n6 = data1[,5]
n7 = data1[,14]
n8 = (data1[,8]+data1[,11]+data1[,12]+data1[,13])
n9 = data1[,15]     


optim5 = optim(y,GetLikSymm,control=list(maxit=1500))
lik5 = optim5$value
conv5 = optim5$convergence
par5 = optim5$par

# Tree 6: (1,(4,(2,3)))

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim6 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik6 = optim6$value
conv6 = optim6$convergence           
par6 = optim6$par

# Tree 7: (2,(1,(3,4))

data1 = matrix(scan("temp4.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim7 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik7 = optim7$value
conv7 = optim7$convergence          
par7 = optim7$par

# Tree 8: (2,(3,(1,4))

data1 = matrix(scan("temp5.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim8 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik8 = optim8$value
conv8 = optim8$convergence               
par8 = optim8$par

# Tree 9: (2,(4,(1,3))

data1 = matrix(scan("temp6.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim9 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik9 = optim9$value
conv9 = optim9$convergence     
par9 = optim9$par

# Tree 10: (3,(1,(2,4))

data1 = matrix(scan("temp7.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim10 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik10 = optim10$value
conv10 = optim10$convergence         
par10 = optim10$par

# Tree 11: (3,(2,(1,4))

data1 = matrix(scan("temp8.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim11 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik11 = optim11$value
conv11 = optim11$convergence       
par11 = optim11$par

# Tree 12: (3,(4,(1,2))

data1 = matrix(scan("temp9.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim12 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik12 = optim12$value
conv12 = optim12$convergence              
par12 = optim12$par

# Tree 13: (4,(1,(2,3))

data1 = matrix(scan("temp10.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim13 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik13 = optim13$value
conv13 = optim13$convergence              
par13 = optim13$par

# Tree 14: (4,(2,(1,3))

data1 = matrix(scan("temp11.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim14 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik14 = optim14$value
conv14 = optim14$convergence              
par14 = optim14$par

# Tree 15: (4,(3,(1,2))

data1 = matrix(scan("temp12.ssa"),ncol=15,byrow=T)

n1 = data1[,1]
n2 = data1[,2]+data1[,3]
n3 = data1[,6]
n4 = data1[,10]
n5 = data1[,4]
n6 = data1[,7]+data1[,9]
n7 = data1[,5]
n8 = data1[,14]
n9 = data1[,8]+data1[,12]
n10 = data1[,11]+data1[,13]
n11 = data1[,15]       

optim15 = optim(y,GetLikAsymm,control=list(maxit=1500))
lik15 = optim15$value
conv15 = optim15$convergence              
par15 = optim15$par

# Now find min of -loglik and check convergence

liks = c(lik1,lik2,lik3,lik4,lik5,lik6,lik7,lik8,lik9,lik10,lik11,lik12,lik13,lik14,lik15)
print(liks)
pars = matrix(c(par1,par2,par3,par4,par5,par6,par7,par8,par9,par10,par11,par12,par13,par14,par15),ncol=4,byrow=T)
print(pars)
mymin = which.min(liks)
#if (mymin==5 || mymin==6 || mymin==8 || mymin==11 || mymin==13) mymin.unroot=3
#else if (mymin==3 || mymin==4 || mymin==9 || mymin==10 || mymin==14) mymin.unroot=2
#else mymin.unroot=1
allconv = conv1+conv2+conv3+conv4+conv5+conv6+conv7+conv8+conv9+conv10+conv11+conv12+conv13+conv14+conv15
w = c(mymin,allconv)
write(w,file="temp.R")

