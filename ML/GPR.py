import math as m
import numpy as np
import torch
import gpytorch
from matplotlib import pyplot as plt

training_iterations = 1000
CUT = 0.5

LABEL = []
yDat = []

# Read data file
TRAINFILE = 'train.dat'
DESIGNFILE = 'Designs.out'
PREDFILE = 'preds.out'

F1 = open(TRAINFILE,'r')
for line in F1:
	spl = line.split()
	if spl[0][0] == '#':
		continue
	else:
		LABEL.append(spl[0])
		yDat.append(float(spl[1]))

F1.close()	


# Read all designs

xTest = []
LABELT = []
F2 = open(DESIGNFILE,'r')
for line in F2:
	spl = line.split()
	xTest.append([float(spl[1]),float(spl[2]),float(spl[3]),float(spl[4])])
	LABELT.append(spl[0])
F2.close()

xDat = []

for i in range(0,len(LABEL)):
	for j in range(0,len(LABELT)):
		if LABEL[i] == LABELT[j]:
			xDat.append(xTest[j])


#===============================================================================#
# Training, heavily based on tutorial from GPyTorch here:https://docs.gpytorch.ai/en/stable/examples/01_Exact_GPs/Simple_GP_Regression.html
train_x = torch.tensor(xDat)
train_y = torch.tensor(yDat)
#yTrain1 = (yTrain1 - min(yDat1))/(max(yDat1)-min(yDat1))

class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        #self.mean_module = gpytorch.means.ZeroMean()
        self.mean_module = gpytorch.means.ConstantMean()
#        self.covar_module = gpytorch.kernels.ScaleKernel(
#            gpytorch.kernels.MaternKernel()
#        )
		
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel()
        )		


    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


likelihood = gpytorch.likelihoods.GaussianLikelihood()
model = ExactGPModel(train_x, train_y, likelihood)

# Find optimal model hyperparameters
model.train()
likelihood.train()

# Use the adam optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=0.1)  # Includes GaussianLikelihood parameters


# "Loss" for GPs - the marginal log likelihood
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
L = []
for i in range(training_iterations):
    optimizer.zero_grad()
    output = model(train_x)
    loss = -mll(output, train_y)
    loss.backward()
    L.append(loss.item())
    print('Iter %d/%d - Loss: %.3f' % (i + 1, training_iterations, loss.item()))
    optimizer.step()

# Set into eval mode
model.eval()
likelihood.eval()


# Make predictions 
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    test_x = torch.tensor(xTest)
    predictions = likelihood(model(test_x))
    mean = predictions.mean
    lower, upper = predictions.confidence_region()

plt.figure()
plt.plot(L)
plt.show()




ER = []


for i in range(0,len(lower)):
	ER.append(upper[i] - lower[i])

	
plt.figure()
plt.hist(ER)
plt.show()

plt.figure()
plt.plot(mean,'*')
plt.show()





LT = np.array(LABELT)
ERP = np.array(ER)


ind = ERP.argsort()


sLT = LT[ind]

sERP = ERP[ind]

NOUT = 5
print('The top errors for the first criteria are:')
for i in range(len(sLT)-1,len(sLT) - NOUT - 1,-1):
	print(str(sLT[i]) + ' with error ' + str(sERP[i]))

print('All candidates below cutoff are:')
for i in range(0,len(LT)):
	if mean[i] < CUT:
		print(str(LT[i]) + ' ' + str(mean[i]))
print('Smallest 20 candidates are:')
min = []
for i in range(0,len(LT)):
	min.append([mean[i],LT[i]])
min.sort()
for i in range(0,20):
	print(str(min[i][1]) + ' ' + str(min[i][0]))

F2 = open(PREDFILE,'w')
for i in range(0,len(LT)):
	F2.write(str(LT[i]) + ' ' + str(float(mean[i])) + ' ' + str(float((upper[i] - lower[i])/2))+ '\n') 
F2.close()









