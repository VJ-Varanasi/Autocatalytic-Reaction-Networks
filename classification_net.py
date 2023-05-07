import numpy as np
import pandas as pd
import sys

from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense

data = pd.read_csv("Data/Kauffman_Data2.csv")


vals = data["RAF"]
ii = np.where(vals == "RAF")[0]
data= data.drop(ii)
data = data.astype(float)
data["f"] = data["p"] * data["edges"]
data = (data - data.min())/(data.max() - data.min())
data =data.dropna(axis = 1)



X_train, X_test, y_train, y_test = train_test_split(data.drop('RAF',axis=1), 
                                                    data['RAF'], test_size=0.30)

print("Data Size: {}".format(len(data)))
print("RAF Sets: {}".format(data["RAF"].sum()))
print("")


l1 = int(sys.argv[1])
l2 = int(sys.argv[2])
l3 = int(sys.argv[3])

dim = len(X_train.columns)

classifier=Sequential()
input_dim = len(data.columns)
classifier.add(Dense(units=l1, input_dim=dim, kernel_initializer='uniform', activation='relu'))
classifier.add(Dense(units=l2, kernel_initializer='uniform', activation='relu'))
classifier.add(Dense(units=l3, kernel_initializer='uniform', activation='relu'))
classifier.add(Dense(units=1, kernel_initializer='uniform', activation='sigmoid'))
classifier.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

BinaryClassifier_Model=classifier.fit(X_train,y_train, batch_size=100 , epochs=10, verbose=1)

pred_y = [1 if x > 0.5 else 0 for x in classifier.predict(X_test)]

true_pos = 0
false_pos = 0
true_neg  = 0
false_neg = 0

y_test = list(y_test)

for i in range(len(y_test)):
  
    if int(y_test[i]) ==1:
        if pred_y[i] == 1:
            true_pos += 1
        else:
            false_neg += 1
    else:
        if pred_y[i] == 1:
            false_pos += 1
        else:
            true_neg += 1


print("True Positive: {}".format(true_pos/len(y_test)))
print("False Positive: {}".format(false_pos/len(y_test)))
print("True Negative: {}".format(true_neg/len(y_test)))
print("False Negative: {}".format(false_neg/len(y_test)))