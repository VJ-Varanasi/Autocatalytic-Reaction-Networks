import numpy as np
import pandas as pd
import sys
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import GridSearchCV


data = pd.read_csv("Data/Kauffman_Data2.csv")

vals = data["RAF"]
ii = np.where(vals == "RAF")[0]
data= data.drop(ii)
data = data.astype(float)


X_train, X_test, y_train, y_test = train_test_split(data.drop('RAF',axis=1), 
                                                    data['RAF'], test_size=0.30)


k = 50
parameters = {"n_neighbors": range(1, k)}
gridsearch = GridSearchCV(KNeighborsRegressor(), parameters)
gridsearch.fit(X_train, y_train)

print(gridsearch.best_estimator_)

train_score= gridsearch.score(X_train, y_train)
print("Train RMSE: {}".format(train_score))

test_score = gridsearch.score(X_test, y_test)
print("Test RMSE: {}".format(test_score))



