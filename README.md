# tRFTargets_Prediction

## 思路
+ 1.利用靶标数据训练模型，并且编写计算tRF与其靶标序列特征的函数  

+ 2.获取全部mRNA序列列表，对于输入的tRF序列进行遍历匹配计算每一个trf-mrna的特征向量  

+ 3.将计算好的特征向量代入训练好的模型中进行预测打分，输出得分高的序列 

+ 4.模型：逻辑回归模型，决策树模型，随机森林模型



## 特征计算（目前得到的特征不太理想，模型学不到太多东西）

+ 1.获取trf-mrna结合位点
+ 2.计算结合能
+ 3.计算AU含量
+ 4.计算结合区结合位点数
+ 5.3'区的结合位点数


## reference

+ [tRForest](https://www.trforest.com/)
+ [文献：tRForest](https://academic.oup.com/nargab/article/4/2/lqac037/6594978)




