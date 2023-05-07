import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

fs = [0.5, 1.5]
data = pd.read_csv("Data/Kauffman_Data_f{}_{}.csv".format(fs[0], fs[1]))


print("fs: {} - {}".format(fs[0], fs[1]))

vals = data["RAF"]
ii = np.where(vals == "RAF")[0]
data= data.drop(ii)
data = data.astype(float)
data["f"] = data["p"] * (data["nodes"] - data["Food Size"])
data["food to node ratio"] =  data["Food Size"]/ data["nodes"]
data= data.drop(columns=['n', 'p'])
data = (data - data.min())/(data.max() - data.min())
data =data.dropna(axis = 1)



X_train, X_test, y_train, y_test = train_test_split(data.drop('RAF',axis=1), 
                                                    data['RAF'], test_size=0.30)

print("Data Size: {}".format(len(data)))
print("RAF Sets: {}".format(data["RAF"].sum()))
print("")


logmodel = LogisticRegression(max_iter = 100000)
logmodel.fit(X_train,y_train)
score = logmodel.score(X_test, y_test)
print("Score: {}".format(score))
print("------------------")

coefs = logmodel.coef_[0]
for i in range(len(coefs)):
    print("{}: {}".format(X_train.columns[i],coefs[i]))
intercept = logmodel.intercept_[0]
print("intercept: {}".format(intercept))


