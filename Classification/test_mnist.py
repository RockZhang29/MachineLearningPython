from sklearn.datasets import fetch_openml
mnist = fetch_openml('mnist_784',version=1)
x,y = mnist["data"],mnist["target"]
# print(x.shape)
# print(y.shape)
x_train,x_test,y_train,y_test = x[:60000],x[60000:],y[:60000],y[60000:]

## 创建二元分类器
from sklearn.linear_model import SGDClassifier
y_train_5 = (y_train == 5)
y_test_5 = (y_test == 5)
sgd_clf = SGDClassifier()
sgd_clf.fit(x_train , y_train_5)

import matplotlib.pyplot as plt
for i in range(0,10):
    some_digit = x[i]
    # some_digit_image = some_digit.reshape(28,28)
    # plt.imshow(some_digit_image,cmap="binary")
    # plt.show()
    ans = sgd_clf.predict(some_digit)
    print(ans)