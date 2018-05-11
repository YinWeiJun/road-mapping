import sys
import json
import math
import numpy
sys.path.append("./")
# from FormateData import FormateData
import xgboost as xgb
from matplotlib import pyplot

# wayTxt真值表
class MappingDatas(object):

    def __init__(self):
        self.testData = []
        self.trainData = []

        self.testLabel = []
        self.trainLabel = []

        self.testRate = []

        self.testId1s = []
        self.testId2s = []

        self.trainId1s = []
        self.trainId2s = []

        self.n_testR = 0
        self.n_testE = 0
        self.testAccuracy = 0

        self.n_testLabel0 = 0
        self.n_testLabel1 = 0

        self.n_trainLabel0 = 0
        self.n_trainLabel1 = 0

        self.n_testErrorData = 0
        self.n_trainErrorData = 0

        self.n_testRPositive = 0
        self.n_testRNegative = 0
        self.N_testRPositive = 0
        self.N_testRNegative = 0

    # 从txt中提取训练数据
    def getTrainDataFromTxt(self, trainTxtPath):
        trainFile = open(trainTxtPath)
        trainLines = trainFile.readlines()

        # wayFile = open(wayTxtPath, 'r')
        # wayLines = wayFile.readlines()
        #
        # osm_id1s = []
        # osm_id2s = []
        #
        # for line in wayLines:
        #     ids = line.split('\t')
        #     osm_id1s.append(ids[1])
        #     osm_id2s.append(ids[0])

        for line in trainLines:
            dataLine = line.split(' ')[-1]
            datas = dataLine.split(',')
            self.trainId2s.append(datas[0])
            self.trainId1s.append(datas[1])
            dataList = []
            dataList.append(float(datas[2])/float(datas[3]))
            dataList.append(float(datas[2])/float(datas[4]))
            dataList.append(float(datas[6]))
            dataList.append(float(datas[5]))
            # div = 1
            # if 1 < float(datas[2]):
            #     div = float(datas[2])
            # dataList.append(float(datas[5])/div)
            self.trainData.append(dataList)

            # if datas[0] in osm_id2s:
            #     index = osm_id2s.index(datas[0])
            #     if datas[1] in osm_id1s[index]:
            #         self.trainLabel.append(1)
            #     else:
            #         self.trainLabel.append(0)

    # 从txt中提取测试数据
    def getTestDataFromTxt(self, testTxtPath):
        testFile = open(testTxtPath)
        testLines = testFile.readlines()

        for line in testLines:
            jsonData = json.loads(line)
            if '0.1' == jsonData['confidence'] and 'way' == jsonData['osm_type']:
                self.testId1s.append(jsonData['osm_id1'])
                self.testId2s.append(jsonData['osm_id2'])
                dataList = []
                dataList.append(float(jsonData['cover1']))
                dataList.append(float(jsonData['cover2']))
                dataList.append(float(jsonData['inc_angle']))
                dataList.append(float(jsonData['test_avgdis']))
                self.testData.append(dataList)

    # 从txt中提取置信度低于confidence的记录
    def getTestDataFromTxtByConfidence(self, confidence, testTxtPath):
        testFile = open(testTxtPath)
        testLines = testFile.readlines()

        for line in testLines:
            jsonData = json.loads(line)
            if confidence > jsonData['confidence'] and 'way' == jsonData['osm_type']:
                self.testId1s.append(jsonData['osm_id1'])
                self.testId2s.append(jsonData['osm_id2'])
                dataList = []
                dataList.append(float(jsonData['cover1']))
                dataList.append(float(jsonData['cover2']))
                dataList.append(float(jsonData['inc_angle']))
                dataList.append(float(jsonData['test_avgdis']))
                self.testData.append(dataList)

    # 丰富的属性，通过didi_sw id联系分别获取didi和sw文件中的属性
    def addAttrToTestDataFromTxt(self, didiTxtPath, swTxtPath):
        if 0 >= len(self.testData):
            return
        n = 0
        didiFile = open(didiTxtPath, 'r')
        swFile = open(swTxtPath, 'r')

        didiLines = didiFile.readlines()
        swLines = swFile.readlines()

        didi_ids = []
        sw_ids = []

        didi_attr = []
        sw_attr = []

        for line in didiLines:
            jsonData = json.loads(line)
            if 'w' not in jsonData['osm_id']:
                continue
            didi_ids.append(jsonData['osm_id'])
            dataList = []
            dataList.append(jsonData['tags']['FuncClass'])
            dataList.append(jsonData['tags']['SpeedClass'])
            dataList.append(jsonData['tags']['kindclass'])
            dataList.append(jsonData['tags']['LaneNum'])
            if 'dual-Carriageway' in jsonData['tags'].keys():
                dataList.append(jsonData['tags']['dual-Carriageway'])#可能为空
            else:
                dataList.append(None)
            if 'structure-link' in jsonData['tags'].keys():
                dataList.append(jsonData['tags']['structure-link'])#可能为空
            else:
                dataList.append(None)
            didi_attr.append(dataList)
        for line in swLines:
            jsonData = json.loads(line)
            if 'w' not in jsonData['osm_id']:
                continue
            sw_ids.append(jsonData['osm_id'])
            dataList = []
            dataList.append(jsonData['tags']['FuncClass'])
            dataList.append(jsonData['tags']['SpeedClass'])
            dataList.append(jsonData['tags']['kindclass'])
            dataList.append(jsonData['tags']['LaneNum'])
            if 'dual-Carriageway' in jsonData['tags'].keys():
                dataList.append(jsonData['tags']['dual-Carriageway'])#可能为空
            else:
                dataList.append(None)
            if 'structure-link' in jsonData['tags'].keys():
                dataList.append(jsonData['tags']['structure-link'])#可能为空
            else:
                dataList.append(None)
            sw_attr.append(dataList)
        for i in range(len(self.testId1s)):
            sw_id = self.testId1s[i]
            didi_id = self.testId2s[i]
            if sw_id not in sw_ids or didi_id not in didi_ids:
                continue
            n += 1
            sw_index = sw_ids.index(sw_id)
            didi_index = didi_ids.index(didi_id)

            dataList = []
            for j in range(len(sw_attr[sw_index])):
                if 0 == j or 1 == j:
                    dataList.append(math.fabs(int(sw_attr[sw_index][j]) - int(didi_attr[didi_index][j])))
                else:
                    if sw_attr[sw_index][j] == None or didi_attr[didi_index][j] == None:
                        dataList.append(numpy.nan)
                    elif sw_attr[sw_index][j] == didi_attr[didi_index][j]:
                        dataList.append(1)
                    else:
                        dataList.append(0)

            self.testData[i] += dataList
        return n

    # 导出训练数据
    def saveTrainData(self, saveBuffer):
        if 0 >= len(self.trainData) or 0>= len(self.trainLabel):
            return
        dtrain = xgb.DMatrix(self.trainData, self.trainLabel)
        dtrain.save_binary(saveBuffer)

    # 导出测试数据
    def saveTestData(self, saveBuffer):
        if 0 >= len(self.testData):
            return
        dtrain = xgb.DMatrix(self.testData)
        dtrain.save_binary(saveBuffer)

    # 加载训练数据
    def loadTrainDataFromBuffer(self, buffer):
        self.dtrain = xgb.DMatrix(buffer)

    # 加载测试数据
    def loadTestDataFromBuffer(self, buffer):
        self.dtest = xgb.DMatrix(buffer)

    # 获取训练模型
    def trainModel(self):
        if len(self.trainLabel) <= 0 or len(self.trainData) <= 0:
            return
        # 加载训练数据
        dtrain = xgb.DMatrix(self.trainData, self.trainLabel)
        # 训练参数
        params = {'max_depth': 3, 'eta': 1, 'silent': 0, 'objective': 'binary:logistic'}
        # 迭代次数
        num_round = 7
        # 训练模型
        bst = xgb.train(params, dtrain, num_round)
        # bst.dump_model('/Users/didi/Documents/confidence01Model.raw.txt')
        self.model = bst
        # bst.save_model('/Users/didi/Documents/confidence01Model.model')
        print('-----------trainModel-------------')

    # 训练模型
    def trainModelByParams(self, params, num_round):
        if len(self.trainLabel) <= 0 or len(self.trainData) <= 0:
            return
        # 加载训练数据
        dtrain = xgb.DMatrix(self.trainData, self.trainLabel)
        # 训练模型
        params = json.loads(params)
        bst = xgb.train(params, dtrain, num_round)
        self.model = bst
        # bst.save_model('/Users/didi/Documents/confidence01Model.model')
        print('-----------getTrainModel-------------')

    # 保存model
    def saveModel(self, modelPath):
        if self.model:
            self.model.save_model(modelPath)

    # 加载模型
    def getTrainModel(self, modelPath):
        bst = xgb.Booster(model_file=modelPath)
        # bst = xgb.Booster()  # init model
        # bst.load_model(modelPath)
        self.model = bst

    # 预测输入的数据
    def predictUserData(self, testData):
        if 0 >= len(testData):
            return 0
        if self.model:
            tData = xgb.DMatrix(testData)
            tLabel = self.model.predict(tData)
            return tLabel

    # 用导入的模型预测数据
    def predictData(self):
        if not self.model:
            return
        dtest_train = xgb.DMatrix(self.testData)
        test_preds = self.model.predict(dtest_train)
        self.testRate = test_preds
        xgb.plot_tree(self.model, num_trees=0, rankdir='LR')
        xgb.plot_importance(self.model)
        xgb.to_graphviz(self.model)
        pyplot.show()

    # 训练模型并预测test数据
    def predictTestData(self):
        if len(self.trainLabel) <= 0 or len(self.testData) <= 0:
            return
        # 加载训练数据
        dtrain = xgb.DMatrix(self.trainData, self.trainLabel)
        # 训练参数
        params = {'max_depth': 3, 'eta': 1, 'silent': 0, 'objective': 'binary:logistic'}
        # params = {
        #     'booster': 'gbtree',
        #     'objective': 'binary:logistic',
        #     'gamma': 0.1,
        #     'max_depth': 4,
        #     'lambda': 3,
        #     'subsample': 0.7,
        #     'colsample_bytree': 0.7,
        #     'min_child_weight': 3,
        #     'silent': 1,
        #     'eta': 1,
        #     'nthread': 4,
        # }
        # 迭代次数
        num_round = 10
        # 训练模型
        bst = xgb.train(params, dtrain, num_round)
        bst.save_model('/Users/didi/Documents/confidence01Model.model')
        # 预测数据
        dtest_train = xgb.DMatrix(self.trainData)
        train_preds = bst.predict(dtest_train)
        # train_predictions = [round(value) for value in train_preds]
        train_predictions = []
        n_predict_r = 0
        n_predict_all = 0
        for value in train_preds:
            if value > 0.98:
                n_predict_all += 1
                train_predictions.append(1)
            else:
                train_predictions.append(0)
        cnt1 = 0
        cnt2 = 0
        for i in range(len(train_predictions)):
            if self.trainLabel[i] == train_predictions[i]:
                cnt1 += 1
                if self.trainLabel[i] == 1:
                    n_predict_r += 1
            else:
                cnt2 += 1

        print("Accuracy: %.2f %% " % (100 * n_predict_r / n_predict_all))
        print("Accuracy: %.2f %% " % (100 * cnt1 / (cnt1 + cnt2)))

        dtest = xgb.DMatrix(self.testData)
        test_preds = bst.predict(dtest)

        # 将概率转换为0／1模型
        # test_predictions = []
        # for value in test_preds:
        #     if 0.9 <= value:
        #         test_predictions.append(1)
        #     else:
        #         test_predictions.append(0)

        self.testRate = test_preds

        # xgb.plot_tree(bst, num_trees=0, rankdir='LR')
        # xgb.plot_importance(bst)
        # xgb.to_graphviz(bst)
        # pyplot.show()
        # pyplot.savefig('foo.png')

    # 将预测结果分类 rateList=[a,b,c], 分4类，0-a a-b b-c c-1
    def classification(self, rateList):
        if len(rateList) <= 0:
            return

        testLabels = []
        for value in self.testRate:
            for i in range(len(rateList) + 1):
                if i == 0 and value < rateList[i]:
                    testLabels.append(i)
                elif i == len(rateList) and rateList[i-1] <= value < 1:
                    testLabels.append(i)
                else:
                    if rateList[i - 1] <= value < rateList[i]:
                        testLabels.append(i)
        self.testLabel = testLabels

    # 获取测试数据正负样本数量
    def getTestDataAccuracy(self, wayTxtPath):
        if len(self.testId1s) <= 0:
            return

        wayFile = open(wayTxtPath, 'r')
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []

        for line in wayLines:
            ids=line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.testId2s)):
            osm_id2 = self.testId2s[i]

            if osm_id2 in osm_id2s:
                index = osm_id2s.index(osm_id2)

                if self.testId1s[i] in osm_id1s[index]:
                    self.n_testLabel1 += 1
                else:
                    self.n_testLabel0 += 1
            else:
                self.n_testErrorData += 1

    # 获取训练样本数据正负样本数量
    def getTrainDataAccuracy(self, wayTxtPath):
        if len(self.trainId1s) <= 0 or len(self.trainId2s) <= 0:
            return

        wayFile = open(wayTxtPath, 'r')
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []

        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.trainId2s)):
            osm_id2 = self.trainId2s[i]

            if osm_id2 in osm_id2s:
                index = osm_id2s.index(osm_id2)

                if self.trainId1s[i] in osm_id1s[index]:
                    self.n_trainLabel1 += 1
                else:
                    self.n_trainLabel0 += 1
            else:
                self.n_trainErrorData += 1

    # 获取预测结果的正确率
    def getTestAnsAccuracy(self, wayTxtPath):
        if len(self.testLabel) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.testId2s)):
            osm_id2 = self.testId2s[i]
            if osm_id2 in osm_id2s:
                index = osm_id2s.index(osm_id2)
                if self.testId1s[i] in osm_id1s[index] and 1 == self.testLabel[i]:
                    self.n_testR += 1
                elif self.testId1s[i] not in osm_id1s[index] and 0 == self.testLabel[i]:
                    self.n_testR += 1
                else:
                    self.n_testE += 1

        self.testAccuracy = self.n_testR / (self.n_testE + self.n_testR)


    # 获取预测结果中正结果站+样本的比例，错误结果占-样本的比例
    def getTestAnsAccuracyEach(self, wayTxtPath):
        if len(self.testLabel) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.testLabel)):
            osm_id2 = self.testId2s[i]
            if 1 == self.testLabel[i]:
                self.N_testRPositive += 1
                if osm_id2 in osm_id2s:
                    index = osm_id2s.index(osm_id2)
                    if self.testId1s[i] in osm_id1s[index]:
                        self.n_testRPositive += 1
            else:
                self.N_testRNegative += 1
                if osm_id2 in osm_id2s:
                    index = osm_id2s.index(osm_id2)
                    if self.testId1s[i] not in osm_id1s[index]:
                        self.n_testRNegative += 1

    # 获取某类预测结果的正确率
    def getTestAnsAccuracyKind(self, kind, wayTxtPath):
        if len(self.testLabel) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        n_R, n_All = 0, 0
        for i in range(len(self.testLabel)):
            osm_id2 = self.testId2s[i]
            if kind == self.testLabel[i]:
                if osm_id2 in osm_id2s:
                    n_All += 1
                    index = osm_id2s.index(osm_id2)
                    if self.testId1s[i] in osm_id1s[index]:
                        n_R += 1
        return (1.0 * n_R) / n_All

    # 0.1抽样训练使用
    def getTrainDataLabel(self, wayTxtPath):
        if len(self.trainData) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.trainId2s)):
            osm_id2 = self.trainId2s[i]
            if osm_id2 in osm_id2s:
                index = osm_id2s.index(osm_id2)
                if self.trainId1s[i] in osm_id1s[index]:
                    self.trainLabel.append(1)
                else:
                    self.trainLabel.append(0)
            else:
                self.trainLabel.append(0)
        print('-------OK---------')

    # 获取测试数据的类别
    def getTestDataLabel(self, wayTxtPath):
        if len(self.trainData) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        for i in range(len(self.testId2s)):
            osm_id2 = self.testId2s[i]
            if osm_id2 in osm_id2s:
                index = osm_id2s.index(osm_id2)
                if self.testId1s[i] in osm_id1s[index]:
                    self.testLabel.append(1)
                else:
                    self.testLabel.append(0)

    # 获取概率大于rate的错误预测项
    def getErrorIds(self, rate, wayTxtPath):
        if len(self.testRate) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        ErrorIds = []
        n_R, n_All = 0, 0
        for i in range(len(self.testRate)):
            rateLabel = self.testRate[i]
            if rateLabel >= rate:
                osm_id2 = self.testId2s[i]
                if osm_id2 in osm_id2s:
                    n_All += 1
                    index = osm_id2s.index(osm_id2)
                    if self.testId1s[i] not in osm_id1s[index]:
                        n_R += 1
                        ErrorIds.append(osm_id2 + ' ' + self.testId1s[i]+'  '+str(self.testRate[i])+'   '+str(self.testData[i]))
        print('---error----',n_R,'---num---',n_All,'-----r_num-----',str(n_All*1.0/len(self.testData))+'----------')
        return ErrorIds

    # 获取概率在rate1，rate2之间的数据的错误结果，kind==1，判断为真的错误项，kind==0，判断为假的错误项
    def getErrorIds2(self, rate1, rate2, wayTxtPath, kind):
        if len(self.testRate) <= 0:
            return

        wayFile = open(wayTxtPath)
        wayLines = wayFile.readlines()

        osm_id1s = []
        osm_id2s = []
        for line in wayLines:
            ids = line.split('\t')
            osm_id1s.append(ids[1])
            osm_id2s.append(ids[0])

        ErrorIds = []
        n_R, n_All = 0, 0
        for i in range(len(self.testRate)):
            rateLabel = self.testRate[i]
            if rate1 < rateLabel <= rate2:
                osm_id2 = self.testId2s[i]
                if osm_id2 in osm_id2s:
                    n_All += 1
                    index = osm_id2s.index(osm_id2)
                    if 1 == kind:
                        if self.testId1s[i] not in osm_id1s[index]:
                            n_R += 1
                            ErrorIds.append(osm_id2 + ' ' + self.testId1s[i]+'  '+str(self.testRate[i])+'   '+str(self.testData[i]))
                    elif 0 == kind:
                        if self.testId1s[i] in osm_id1s[index]:
                            n_R += 1
                            ErrorIds.append(osm_id2 + ' ' + self.testId1s[i]+'  '+str(self.testRate[i])+'   '+str(self.testData[i]))
        print('---error----',n_R,'---num---',n_All,'-----r_num-----',str(n_All*1.0/len(self.testData))+'----------')
        return ErrorIds
