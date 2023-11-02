#!/usr/bin/env python3
# coding: utf-8
# author: zhang.siwen
# references:
#     https://scikit-learn.org/stable/index.html
#     https://chrisalbon.com/#code_machine_learning
#     https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
#     https://chrisalbon.com/code/machine_learning/model_selection/find_best_preprocessing_steps_during_model_selection/


# common libs
import argparse
import sys
import os
import pandas as pd
import numpy as np
import joblib
import warnings
warnings.filterwarnings('ignore')

# scale and decomposition libs
from sklearn.preprocessing import StandardScaler, MinMaxScaler, PowerTransformer, RobustScaler, QuantileTransformer, label_binarize
from sklearn.feature_selection import SelectKBest
from sklearn.decomposition import PCA

# classifier libs
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier

# model training
from sklearn import metrics
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline, FeatureUnion

# plot lib
import matplotlib.pyplot as plt
import seaborn as sns
# plt.style.use('seaborn-pastel')


# fixed parameters
FEATURE_NUMBER = 500



def GetArgs():
    parser = argparse.ArgumentParser(description='Find best classify method and parameters', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--matrix', help='matrix to be classified, with rows as features and columns as samples', action='store', dest='matrix', required=True)
    required.add_argument('--stats', help='for knowning answer, output confusion matrix statistics of all classifier,\n  for not knowning answer, output prediction result', action='store', dest='stats', required=True)
    # model training and testing 
    fitmodel = parser.add_argument_group('Model training and testing arguments')
    fitmodel.add_argument('--groups', help='sample groups in TSV format, \n  first column is sample name, second column is sample group', default=None, action='store', dest='groups', required=False)
    fitmodel.add_argument('--model', help='path to read model for testing matrix, DEFAULT "%(default)s"', default=None, action='store', dest='model', required=False)
    fitmodel.add_argument('--save', help='path to save model for re-use, DEFAULT "%(default)s"', default=None, action='store', dest='save', required=False)
    fitmodel.add_argument('--train', help='train set ratio, DEFAULT "%(default)s"', default=1, type=float, action='store', dest='train', required=False)
    fitmodel.add_argument('--method', help='classifier method, if multiple method used, use comma "," to seperate, \n  DEFAULT "%(default)s" \n  types: {}'.format(str(PredictorClasses.keys())), default='rf', choices=list(PredictorClasses.keys())+['auto'], action='store', dest='method', required=False)
    scoring_method = ['none', 'roc_auc','roc_auc_*','accuracy','f1','f1_*','recall','precision']
    fitmodel.add_argument('--scoring', help='scoring method for assessing parameters in GridSearchCV, \n  use "none" for build model with default parameter, \n  DEFAULT "%(default)s" \n  types: {}'.format(str(scoring_method)), default='roc_auc', choices=scoring_method, action='store', dest='scoring', required=False)
    # preprocess
    preprocess = parser.add_argument_group('Preprocess arguments')
    preprocess.add_argument('--transform', help='matrix log2 transform, DEFAULT "%(default)s"', default='log2', choices=list(NumberTransform.keys())+['none'], action='store', dest='transform', required=False)
    preprocess.add_argument('--genes', help='gene list to train and test, ignore for only testing situation, DEFAULT "%(default)s"', default=None, action='store', dest='genes', required=False)
    preprocess.add_argument('--refgenes', help='inner reference gene list for calculate relative values, ignore for only testing situation, accept gene name separated by comma "," or a file with each gene in a row, \n  DEFAULT "%(default)s"', default=None, action='store', dest='refgenes', required=False)
    preprocess.add_argument('--corr', help='if sum( corr(Sample,RestSamples)<X ) < Y*TotalSample, remove sample, ignore for only testing situation, \n  correlation method choices: [pearson, kendall, spearman], \n  use "none" for not filtering, \n  format: X,Y,method, separated by comma ",", DEFAULT "%(default)s"', default='0.5,0.1,spearman', action='store', dest='corr', required=False)
    preprocess.add_argument('--feature', help='row-wise sum( feature>=X ) >= Y*samples, keep feature, ignore for only testing situation,\n  correlation method choices: [pearson, kendall, spearman], \n  use "none" for not filtering, \n  format: X,Y,method, separated by comma ",", DEFAULT "%(default)s"', default='5,0.1,spearman', action='store', dest='feature', required=False)
    preprocess.add_argument('--scaler', help='scale data, ignore for only testing situation, DEFAULT "%(default)s" \n  types: {}'.format(str(ScaleClasses.keys())), default='none', choices=list(ScaleClasses.keys())+['none'], action='store', dest='scaler', required=False)
    # other parameters
    other = parser.add_argument_group('Other parameters for model training')
    other.add_argument('--seed', help='seed for model building, should less than 2**32-1,\n  use "-1" for random seed,\n  DEFAULT "%(default)s"', default=1234, type=int, action='store', dest='random_state', required=False)
    # other.add_argument('--seed', help='seed for model building, \n  use "-1" for random seed,\n default %(default)s', default=1234, type=int, choices=range(0,2**32), metavar="[0,2**32-1]", action='store', dest='random_state', required=False)
    other.add_argument('--threads', help='threads for model tuning, \n  use "-1" for all threads available,\n  DEFAULT "%(default)s"', default=1, type=int, action='store', dest='threads', required=False)
    # scale and decomposition
    # decomposition = parser.add_argument_group('Scale and decompositon arguments')
    # decomposition.add_argument('--decom', help='data decomposition (and keep one feature in a set of highly correlated features) by default, DEFAULT "%(default)s" \n types: {}'.format(str(DecompositionClasses)), default='none', choices=list(DecompositionClasses.keys())+['none', 'auto'], action='store', dest='decom', required=False)
    # classify
    
    usage = '''Usage:

For training only:
    %(prog)s --matrix merged.htseq.count --stats output.df --groups sample.list --save output.model [--genes gene.list --refgenes reference_gene.list --corr 0.5,0.1,spearman --feature 5,0.1,spearman --scaler standard --transform log2 ]

For training and testing:
    %(prog)s --matrix merged.htseq.count --stats output.df --groups sample.list --save output.model --train 0.7 [--genes gene.list --corr 0.5,0.1,spearman --feature 5,0.1,spearman --scaler standard --transform log2 ]

For training with default parameter of classifier:
    %(prog)s --matrix merged.htseq.count --stats output.df --groups sample.list --save output.model --scoring none

For testing only and not knowning answer:
    %(prog)s --matrix merged.htseq.count --stats output.df --model trained.model 

For testing only and knowing answer:
    %(prog)s --matrix merged.htseq.count --stats output.df --groups sample.list --model trained.model 

'''
    parser.epilog = usage
    args = parser.parse_args()
    if args.model is None and args.groups is None:
        sys.exit('At least give one of these arguments: --group --model')
    if args.random_state ==  -1:
        args.random_state = None

    return args



NumberTransform = {
    'log2': np.log2,
    'sqrt': np.sqrt,
}

ScaleClasses = {
    'standard': StandardScaler, #might behave badly if the individual features do not more or less look like standard normally distributed data
    'minmax': MinMaxScaler,
    'power': PowerTransformer,
    'robust': RobustScaler,
    'quantile': QuantileTransformer,
}

DecompositionClasses = {
    #'pca': PCA,
    'kbest': SelectKBest,
}

DecompositionArgs = {
    #'pca': {'n_components': np.arange(0, 1, 0.05)},
    'kbest': {'score_func': ['f_classif', 'chi2'],
              'k': np.arange(1, 30, dtype=int)},
}

PredictorClasses = {
    'svm': SVC,
    'logistic': LogisticRegression,
    'knn': KNeighborsClassifier,
    'rf': RandomForestClassifier,
    'dt': DecisionTreeClassifier,
    'adaboost': AdaBoostClassifier,
    'naivebayes': GaussianNB,
    'sgd': SGDClassifier,
    'dnn': MLPClassifier,
}

PredictorArgs = {
    'svm': {'C': np.logspace(-3, 5, 9), 
            'gamma': np.logspace(-9, 3, 13),
            'kernel': ['linear', 'poly', 'rbf', 'sigmoid'] },
    'logistic': {'C': np.arange(0.1, 10, 0.1), 
                 'penalty':['l1','l2'] },
    'knn': {'n_neighbors': np.arange(5, 50, 5, dtype=int), 
            'leaf_size': np.arange(3, 30, dtype=int) },
    'rf': {'max_depth': np.arange(2, 10, dtype=int), 
           # 'max_features': np.arange(10, 100, dtype=int), 
           # 'max_leaf_nodes': np.arange(3, 100, dtype=int),
           'n_estimators': np.arange(10, 100, 10, dtype=int)},
    'dt': {'criterion': ['gini', 'entropy', 'log_loss'], 
           'max_depth': np.arange(5, 50, dtype=int), 
           'max_features': np.arange(5, 100, 5, dtype=int), 
           'max_leaf_nodes': np.arange(3, 30, 3, dtype=int) },
    'adaboost': {'n_estimators': np.arange(5, 100, 5, dtype=int),
                 'learning_rate': np.arange(0.01, 1, 0.01) },
    'naivebayes': {'var_smoothing': np.logspace(0, -9, 20) },
    'sgd': {'loss': ['hinge', 'log_loss', 'perceptron'],
             'penalty': ['l2', 'l1', 'elasticnet'],
             'alpha': np.arange(0.01, 1, 0.01),
             'max_iter': np.arange(100, 1000, 100, dtype=int),
             'learning_rate': ['constant', 'optimal', 'invscaling', 'adaptive'] },
    'dnn': {'hidden_layer_sizes': np.arange(10, 100, 10, dtype=int),
            'activation': ['identity', 'logistic', 'tanh', 'relu'],
            'alpha': np.arange(0.01, 1, 0.01),
            'learning_rate_init': np.arange(0.1, 1, 0.1),
            'max_iter': np.arange(100, 1000, 100, dtype=int) }
}
# parameter range can be found by validation_curve

StepClass = {
    'scale': ScaleClasses,
    'decom': DecompositionClasses,
    # 'classify': PredictorClasses,
}

StepArgs = {
    # 'scale': ScaleArgs,
    'decom': DecompositionArgs,
    'classify': PredictorArgs,
}

# get arguments
# set(PredictorClasses['rf'].__init__.__code__.co_varnames)



# matrix value transform
def MatrixTransform(matrix, method='log2'):
    print("Transform matrix: "+method)
    if method != 'none':
        if method == 'log2':
            matrix = matrix+1 # log2(x+1) for 'log2', sqrt(x) for 'sqrt'
        matrix = NumberTransform[method](matrix)
    return matrix



# if give reference gene list, calculate relative expression levels
def RelativeValuesByReferenceGenes(matrix, refgene_list, calculate_method='median', log2level=True):
    print('Relative value calculation start.')
    print('Number of sample: {}\nNumber of feature: {}'.format(str(matrix.shape[1]),str(matrix.shape[0])))
    refgene_intersect = []
    for i in refgene_list:
        if i in matrix.index:
            refgene_intersect.append(i)
        else:
            print("Reference gene NOT in matrix: "+i)
    # calculate_method choices: ['median', 'mean', 'mad']
    if calculate_method == 'median':
        refgene_standards = matrix.loc[refgene_intersect].median()
    if calculate_method == 'mean':
        refgene_standards = matrix.loc[refgene_intersect].mean()
    if calculate_method == 'mad':
        refgene_standards = matrix.loc[refgene_intersect].mad()

    zero_ref_samples = refgene_standards[refgene_standards==0].index
    if len(zero_ref_samples) != 0:
        sys.stdout("Samples with zero {} values of reference genes:\n{}".format(calculate_method, '\n'.join(zero_ref_samples.to_list())))
    
    if log2level:
        matrix = matrix.drop(columns=zero_ref_samples) - refgene_standards.drop(zero_ref_samples)
    else:
        matrix = matrix.drop(columns=zero_ref_samples) / refgene_standards.drop(zero_ref_samples)
    print('After calculate relative values,')
    print('Number of sample: {}\nNumber of feature: {}'.format(str(matrix.shape[1]),str(matrix.shape[0])))

    return matrix



# preprocess, filter samples
def CheckCorrelation(matrix, corr_score=0.5, corr_ratio=0.1, corr_method='spearman', drop=True):
    # use correlation to get samples that are similar (>0.5) with at least 10% of total samples
    correlation_matrix = matrix.corr(method=corr_method)
    failed_samples = correlation_matrix[ correlation_matrix.ge(corr_score).sum() < corr_ratio*correlation_matrix.shape[0] ].index.to_list()
    if drop == True:
        matrix = matrix.drop(failed_samples, axis=1)
    
    return matrix, correlation_matrix, failed_samples



# preprocess, filter features
def CheckFeatures(matrix, feature_count=5, feature_ratio=0.1, feature_method='spearman', drop=True):
    failed_features = []
    if feature_count != -1:
        failed_features = matrix[ matrix.ge(feature_count).sum(axis=1) >= feature_ratio*matrix.shape[1] ].index.to_list()
    if matrix.shape[0] > FEATURE_NUMBER :
        correlation_matrix = matrix.transpose().corr(method=feature_method)
        upper_matrix = correlation_matrix.where(np.triu(np.ones(correlation_matrix.shape), k=1).astype(np.bool))
        failed_features = failed_features + [column for column in upper_matrix.columns if any(upper_matrix[column] > 0.95)]
    if drop == True and failed_features != list() :
        matrix = matrix.drop(failed_features, axis=0)
    
    return matrix, failed_features



# main function of preprocess data WITH y info
def PreprocessMatrix(matrix, group, key='subtype', 
                     feature_count=5, feature_ratio=1, feature_method='spearman', 
                     corr_score=0.5, corr_ratio=0.02, corr_method='spearman'):
    '''
    matrix: name	s1	s2	s3	s4
            A1BG	0	0	0	0
            A1CF	0	0	0	0
            A2M	0	0	19	1
            A2ML1	0	0	0	0
    group:  sample	key	other_columns_will_not_be_used
            s2	t1	1
            s2	t2	2
            s3	t3	3
            s4	t4	4
    '''
    print('Preprocessing start.')
    print('Number of sample: {}\nNumber of feature: {}'.format(str(matrix.shape[1]),str(matrix.shape[0])))
    try:    
        group = group.loc[matrix.columns.to_list()]
        group.reindex(matrix.columns.to_list())
    except Exception as e: 
        print('Sample not in group: '+ str(e))
        shared_samples = [idx for idx in matrix.columns.to_list() if idx in group.index.to_list()]
        group = group.loc[shared_samples]
        group.reindex(shared_samples)
        matrix = matrix[shared_samples]
    print('Number of sample: {}\nNumber of feature: {}'.format(str(matrix.shape[1]),str(matrix.shape[0])))
    
    if feature_count != None and feature_ratio != None:
        matrix, _ = CheckFeatures(matrix, feature_count=feature_count, feature_ratio=feature_ratio, feature_method=feature_method, drop=True)
        print('Filter features: sum( feature<{} ) > {}*samples'.format(str(feature_count),str(feature_ratio)))
        print('Number of feature: {}'.format(str(matrix.shape[0])))
    if corr_score != None and corr_ratio!=None:
        matrix, _, _ = CheckCorrelation(matrix, corr_score=corr_score, corr_ratio=corr_ratio, corr_method=corr_method, drop=True)
        print('Filter samples: sum( corr(Sample,RestSamples)<{} ) < {}*TotalSample'.format(str(corr_score),str(corr_ratio)))
        print('Number of feature: {}'.format(str(matrix.shape[1])))
    
    xdata = matrix.transpose() #.to_numpy()
    ydata = group[key] #pd.Categorical(group[key]).codes
    ylabel = dict(sorted(zip(pd.Categorical(group[key]).codes, group[key])))
    print('Group transform: '+str(ylabel))
    features = matrix.index.to_list()
    
    return xdata, ydata, ylabel, features



# main function of preprocess data WITHOUT y info
def PreprocessMatrixForTesting(matrix, ylabel, group=None, key='subtype', features=None):
    print('Preprocessing start.')
    print('Number of sample: {}\nNumber of feature: {}'.format(str(matrix.shape[1]),str(matrix.shape[0])))
    if features is not None:
        matrix = matrix[matrix.index.isin(list(features))].reindex(list(features))
        print('Filter features.\nNumber of feature: {}'.format(str(matrix.shape[0])))
    if group is not None:
        # ylabel: {0:'t1', 1:'t2'}
        # label_reverse:  {'t1':0, 't2':1}
        shared_samples = [idx for idx in matrix.columns.to_list() if idx in group.index.to_list()]
        group = group.loc[shared_samples]
        group.reindex(shared_samples)
        matrix = matrix[shared_samples]
        # label_reverse = {v:k for k,v in ylabel.items()}
        # ydata = [label_reverse[i] for i in group[key]] 
        ydata = group[key]
    else:
        ydata = None
    xdata = matrix.transpose() #.to_numpy()
    
    return xdata, ydata



# recognize steps and arguments
def StepRecognize(category, key):
    if key == 'none':
        key = None
    elif key == 'auto':
        key = FeatureUnion([(k,v()) for k,v in StepClass[category].items() ])
    else:
        key = StepClass[category][key]()
    
    return key



def GridParameterFlatten(list_args, key):
    '''
    flatten list_args and add prefix
    list_args = [{'C': array([1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03, 1.e+04, 1.e+05,
       1.e+06, 1.e+07, 1.e+08, 1.e+09, 1.e+10]), 'gamma': array([1.e-09, 1.e-08, 1.e-07, 1.e-06, 1.e-05, 1.e-04, 1.e-03, 1.e-02,
       1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03])}]
    '''
    args = {}
    for k,v in list_args[0].items():
        args[key+'__'+k] = v
    return args



# get parameters of pipelines
def GridPipelineParamaters(pipeline):
    '''
    pipeline = [('scale',
                  FeatureUnion(transformer_list=[('standard', StandardScaler()),
                                                 ('minmax', MinMaxScaler()),
                                                 ('quantile', QuantileTransformer())])),
                 ('decomposition',
                  FeatureUnion(transformer_list=[('pca', PCA()), ('kbest', SelectKBest())])),
                 ('classifier', SVC())]
    '''
    parameters = {}
    for prefix,value in pipeline:
        if prefix == 'scale':
            continue
        if value.__class__ == FeatureUnion:
            for key,method in value.transformer_list:
                if prefix == 'decomposition':
                    args = [DecompositionArgs[k] for k,v in DecompositionClasses.items() if v == method.__class__]
                else:
                    args = [PredictorArgs[k] for k,v in PredictorClasses.items() if v == method.__class__]
                parameters.update(GridParameterFlatten(args, prefix+'__'+key))
        else:
            if prefix == 'decomposition':
                args = [DecompositionArgs[k] for k,v in DecompositionClasses.items() if v == value.__class__]
            else:
                args = [PredictorArgs[k] for k,v in PredictorClasses.items() if v == value.__class__]
            parameters.update(GridParameterFlatten(args, prefix))
    
    return [parameters]
    


# build model
def ModelBuilding(ylabel, method='rf', scaler=None, decom=None, cv=5, scoring='roc_auc', random_state=1234, n_jobs=-1):
    Predictor = PredictorClasses[method]
    if method in ['knn', 'naivebayes']:
        predictor = Predictor()
    # elif method == 'svm':
        # predictor = Predictor(probability=True)
        # https://scikit-learn.org/stable/modules/svm.html#scores-probabilities
        # the probability estimates may be inconsistent with the scores:
        #     - the “argmax” of the scores may not be the argmax of the probabilities
        #     - in binary classification, a sample may be labeled by predict as belonging to the positive class even if the output of predict_proba is less than 0.5; and similarly, it could be labeled as negative even if the output of predict_proba is more than 0.5.
        # https://stackoverflow.com/questions/17017882/scikit-learn-predict-proba-gives-wrong-answers
    else:
        if method=='svm' and len(ylabel)>2:
            predictor = Predictor(probability=True, break_ties=True, decision_function_shape='ovr', random_state=random_state)
        else:
            predictor = Predictor(probability=True, random_state=random_state)
    
    # build pipeline
    pipeline = [(i,j) for i,j in [('scale', scaler), ('decomposition', decom), ('classifier', predictor)] if j!=None]
    
    if scoring == 'none':
        # model = Pipeline(pipeline)
        # set search space to pipeline default parameters
        param_grid = [{k:[v] for k,v in Pipeline(pipeline).get_params().items() if k.startswith('scale') or k.startswith('decomposition') or k.startswith('classifier')}]
        scoring = 'roc_auc'
    else:
        param_grid = GridPipelineParamaters(pipeline)
    
    model = GridSearchCV(Pipeline(pipeline), param_grid=param_grid, cv=cv, scoring=scoring, n_jobs=n_jobs)
    
    return model



# model predict result df
def ModelPredict(fit_models, x, sample_names, y=None):
    # result, {method:[pred1, pred2, ...]}
    model_result = {}
    for method in fit_models.keys():
        model = fit_models[method]
        y_predict = model.predict(x)
        model_result[method] = y_predict #[ylabel[i] for i in y_predict]
    if y is not None:
        model_result['true_label'] = y
    model_predictDF = pd.DataFrame.from_dict(model_result, orient='columns', dtype=str)
    model_predictDF.index = sample_names
    
    return model_predictDF



# calculate confusion matrix and transform to dataframe with TP,FP,FN,TN,sensitivity,specificity,f1score,accuracy calculated
# reference: https://stackoverflow.com/questions/27959467/how-to-find-tp-tn-fp-and-fn-values-from-8x8-confusion-matrix
def FormatConfusionMatrix(y_true, y_pred, ylabel):
    TFcols = ['TP', 'FP', 'FN', 'TN', 'Accuracy', 'Specificity']
    creteria = ['Precision', 'Recall', 'F1score', 'Support']
    
    # confusion matrix
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(y_true, y_pred), columns=ylabel.values())
    confusion_matrix.index = ylabel.values()
    
    # precision report
    report = metrics.classification_report(y_true, y_pred, target_names=ylabel.values(), output_dict=True)
    reportDf = pd.DataFrame(report).transpose()
    reportDf.columns = creteria
    
    # result df
    confusion_matrix = confusion_matrix.join(pd.DataFrame(index=ylabel.values(), columns=TFcols+creteria))
    # confusion_matrix = pd.merge(confusion_matrix, pd.DataFrame(index=ylabel.values(), columns=TFcols), left_index=True, right_index=True)
    confusion_matrix.loc[ylabel.values(),creteria] = reportDf.loc[ylabel.values(),creteria]
    # calculate TP,FP,FN,TN,sensitivity,specificity,f1score,accuracy for each category
    total = confusion_matrix.loc[:, ylabel.values()].sum().sum()
    for key in ylabel.values():
        TP = confusion_matrix.loc[key, key]
        FP = confusion_matrix.loc[ylabel.values(), key].sum() - TP
        FN = confusion_matrix.loc[key, ylabel.values()].sum() - TP
        TN = total - FP - TP - FN
        specificity = TN / (TN+FN)
        accuracy = (TP+TN) / total
        confusion_matrix.loc[key, TFcols] = [TP, FP, FN, TN, accuracy, specificity]
    
    return confusion_matrix



# model fit by tuning best hyper-parameters
def ModelTuningGrid(x_train, y_train, ylabel, sample_names, 
                    scaler='none', decom='none', method=['rf'], scoring='roc_auc',  random_state=1234, n_jobs=-1):
    
    # cv definition
    if len(y_train) < 100:
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state) 
        #StratifiedShuffleSplit(n_splits=5, test_size=0.3, random_state=random_state)
    else:
        cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=random_state) 
        #StratifiedShuffleSplit(n_splits=10, test_size=0.3, random_state=random_state)
    
    # scale and decomposition
    scaler = StepRecognize('scale', scaler)
    decom = StepRecognize('decom', decom)
    
    # fit for each classifier
    train_stats = pd.DataFrame()
    test_stats = pd.DataFrame()
    fit_models = {}
    
    # get method search space
    if method == ['auto']:
        method_grid = list(PredictorClasses.keys())
    else:
        method_grid = method
    print('Training method: {}'.format(','.join([i for i in method_grid])))
    
    # fit models
    for meth in method_grid:
        try:
            search_model = ModelBuilding(ylabel, method=meth, scaler=scaler, decom=decom, cv=cv, scoring=scoring, random_state=random_state, n_jobs=n_jobs)
        except:
            search_model = ModelBuilding(ylabel, method=meth, scaler=scaler, decom=decom, cv=cv, scoring='accuracy', random_state=random_state, n_jobs=n_jobs)
        search_model.fit(x_train, y_train)
        
        # initialize model using best hyper-parameters and fit data
        # if scoring == 'none':
        #     best_param = "default_parameter"
        #     model = search_model
        # else:
        best_param = search_model.best_params_
        model = search_model.estimator.set_params(**best_param)
        model.fit(x_train, y_train)

        # confusion matrix calculation
        train_confusion_matrix = FormatConfusionMatrix(y_train, model.predict(x_train), ylabel)
        try:
            train_confusion_matrix['ROC_AUC'] = round(metrics.roc_auc_score(y_train, model.decision_function(x_train)), 4)
        except:
            train_confusion_matrix['ROC_AUC'] = round(metrics.roc_auc_score(y_train, model.predict_proba(x_train)[:, 1]), 4)
        classifier_name = str(model.steps[-1][-1].__class__).strip("'>").split('.')[-1]
        train_confusion_matrix['Classifier'] = classifier_name
        train_confusion_matrix['Estimator'] = str(best_param)
        train_confusion_matrix.index = [classifier_name+'-'+i for i in train_confusion_matrix.index]
        # save best models
        fit_models[meth] = model
        train_stats = train_stats.append(train_confusion_matrix)
        print('{0}: Train for {0} is done.'.format(meth))
    
    model_train_predictDF = ModelPredict(fit_models=fit_models, x=x_train, sample_names=sample_names, y=y_train)
    print('All training done.')
    
    return fit_models, model_train_predictDF, train_stats



# model test by model
def ModelTestingOnly(fit_models, x_test, ylabel, sample_names, y_test=None):

    # store test stats dataframe
    test_stats = pd.DataFrame()
    # fit models
    for method in fit_models.keys():
        model = fit_models[method]
        if y_test != None:
            # confusion matrix calculation
            y_predict = model.predict(x_test)
            classifier_name = str(model.steps[-1][-1].__class__).strip("'>").split('.')[-1]
            test_confusion_matrix = FormatConfusionMatrix(y_test, y_predict, ylabel)
            try:
                test_confusion_matrix['ROC_AUC'] = round(metrics.roc_auc_score(y_test, model.decision_function(x_test)), 4)
            except:
                test_confusion_matrix['ROC_AUC'] = round(metrics.roc_auc_score(y_test, model.predict_proba(x_test)[:, 1]), 4)
            test_confusion_matrix['Classifier'] = classifier_name
            test_confusion_matrix['Estimator'] = str(fit_models[method].get_params())
            test_confusion_matrix.index = [classifier_name+'-'+i for i in test_confusion_matrix.index]
            test_stats = test_stats.append(test_confusion_matrix)
            print('{0}: Test for {0} is done.'.format(method))
    model_predictDF = ModelPredict(fit_models=fit_models, x=x_test, sample_names=sample_names, y=y_test)
    print('All testing done.')
    
    return model_predictDF, test_stats



# multiclass ROC metrics calculation
def MulticlassROCcalculation(model, x_true, y_true, ylabel):

    # format y values
    y_true = label_binarize(y_true, classes=list(ylabel.keys()))
    y_pred = model.decision_function(x_true)

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in list(ylabel.keys()):
        fpr[i], tpr[i], _ = metrics.roc_curve(y_true[:, i], y_pred[:, i])
        roc_auc[i] = metrics.auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = metrics.roc_curve(y_true.ravel(), y_pred.ravel())
    roc_auc["micro"] = metrics.auc(fpr["micro"], tpr["micro"])

    return fpr, tpr, roc_auc



# plot multiclass ROC
# reference: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
def MulticlassROCPlot(fpr, tpr, roc_auc, ylabel):
    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in list(ylabel.keys())]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in list(ylabel.keys()):
        mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= len(ylabel.keys())

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = metrics.auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    plt.figure()
    for i,color in zip(["micro","macro"], ["deeppink","navy"]):
        plt.plot(
            fpr[i],
            tpr[i],
            label="{0}-average ROC curve (area = {1:0.2f})".format(i, roc_auc[i]),
            color=color,
            linestyle=":",
            linewidth=2,
        )
    
    #colors = cycle(["aqua", "darkorange", "cornflowerblue"])
    colors = sns.color_palette("pastel", len(list(ylabel.keys())))
    for i,color in zip(ylabel.keys(), colors):
        plt.plot(
            fpr[i],
            tpr[i],
            color=color,
            linewidth=2,
            label="ROC curve of class {0} (area = {1:0.2f})".format(ylabel[i], roc_auc[i]),
        )

    plt.plot([0, 1], [0, 1], "k--", lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    #plt.title("Receiver operating characteristic to multiclass")
    plt.legend(loc="lower right")
    #plt.show()
    return plt



# binary ROC plot
def BinaryROCPlot(fit_models, x, y, ylabel):
    
    fig, ax_roc = plt.subplots(1, 1, figsize=(8, 8))

    for name,model in fit_models.items():
        # try:
        #     y_decision = model.decision_function(x)
        # except:
        #     y_decision = model.predict_proba(x)
        # metrics.RocCurveDisplay.from_predictions(y, y_decision, ax=ax_roc, name=name)
        metrics.RocCurveDisplay.from_estimator(model, x, y, ax = ax_roc, name=name)

    ax_roc.set_title("Receiver Operating Characteristic (ROC) curves")
    ax_roc.grid(linestyle="--")

    plt.plot([0, 1], [0, 1], "k--", lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend(loc="lower right")
    
    return plt



# general function for ROC plot
def ROCforModels(fit_models, x, y, ylabel, outpath=None):
    if y == None:
        return
    if len(ylabel) == 2:
        plot = BinaryROCPlot(fit_models, x, y, ylabel)
        if outpath != None:
            plot.savefig(outpath+'.binaryROC.pdf')
    else:
        for method in fit_models.keys():
            fpr, tpr, roc_auc = MulticlassROCcalculation(fit_models[method], x, y, ylabel)
            plot = MulticlassROCPlot(fpr, tpr, roc_auc, ylabel)
            if outpath != None:
                plot.savefig(outpath+'.multiclassROC.'+method+'.pdf')
    return



def Main():
    args = GetArgs()
    count_matrix = pd.read_csv(args.matrix, sep='\t', low_memory=False, index_col=0)
    print('Number of sample: {}\nNumber of feature: {}'.format(str(count_matrix.shape[1]),str(count_matrix.shape[0])))
    stats_path = args.stats.strip('.tsv') if args.stats.endswith('.tsv') else args.stats
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)

    if args.groups is not None:
        sample_group = pd.read_csv(args.groups, sep='\t', low_memory=False, index_col=0, dtype={'sample':str})
        sample_group.index = [str(i) for i in sample_group.index.to_list()]
    else:
        sample_group = None
    
    # for training or for testing
    if args.model is None: # training and/or testing
        # transform matrix
        count_matrix = MatrixTransform(count_matrix, method=args.transform)
        # calculate relative expression level by reference genes
        if args.refgenes is not None:
            if os.path.isfile(args.refgenes):
                refgene_list = [line.strip() for line in open(args.refgenes, 'r')]
            elif ',' in args.refgenes:
                refgene_list = args.refgenes.strip().split(',')
            else:
                refgene_list = args.refgenes.strip()
                if refgene_list in count_matrix.index:
                    refgene_list = [refgene_list]
                else:
                    sys.exit('For reference gene list, provide either a string separated genes with comma "," OR a file with each gene in one row.')
            refgene_cal = 'mean'
            count_matrix = RelativeValuesByReferenceGenes(count_matrix, refgene_list, calculate_method=refgene_cal, log2level=True)
        else:
            refgene_list = None
            refgene_cal = None
        # keep provided genes only
        if args.genes is not None:
            gene_list = pd.read_csv(args.genes, sep='\t', low_memory=False, index_col=0)
            count_matrix = count_matrix[count_matrix.index.isin(gene_list.index.to_list())]
            print('Filter features: ')
            print('Number of sample: {}\nNumber of feature: {}'.format(str(count_matrix.shape[1]),str(count_matrix.shape[0])))
        # preprocessing
        if args.feature.strip() == 'none':
            feature_count, feature_ratio, feature_method = None, None, None
        else:
            feature_count, feature_ratio, feature_method = args.feature.strip().split(',')
        if args.corr.strip() == 'none':
            corr_score, corr_ratio, corr_method = None, None, None
        else:
            corr_score, corr_ratio, corr_method = args.corr.strip().split(',')
        xdata, ydata, ylabel, features = PreprocessMatrix(matrix=count_matrix, group=sample_group, key='Group', 
            feature_count=feature_count, feature_ratio=feature_ratio, feature_method=feature_method,
            corr_score=corr_score, corr_ratio=corr_ratio, corr_method=corr_method)
        
        # parse method search space
        fit_methods = str(args.method).strip().split(',')
        # fit model and summarize stats
        if args.train != 1: # both training and testing
            x_train, x_test, y_train, y_test = train_test_split(xdata, ydata, random_state=args.random_state, train_size=args.train, stratify=ydata)
            fit_models, model_train_predictDF, train_stats = ModelTuningGrid(x_train=x_train, y_train=y_train,
                ylabel=ylabel, sample_names=x_train.index.to_list(),  
                scaler=args.scaler, decom='none', method=fit_methods, scoring=args.scoring, random_state=args.random_state, n_jobs=args.threads)
            model_test_predictDF, test_stats = ModelTestingOnly(fit_models=fit_models, x_test=x_test, ylabel=ylabel, sample_names=x_test.index.to_list(), y_test=y_test.to_list())
            model_test_predictDF.to_csv(stats_path+'.test_predict.tsv', sep='\t', index=True, mode='w')
            test_stats.to_csv(stats_path+'.test_stats.tsv', sep='\t', index=True, mode='w')
            ROCforModels(fit_models=fit_models, x=x_test, y=y_test.to_list(), ylabel=ylabel, outpath=stats_path+'.test_data')
        else: # only training
            x_train, x_test, y_train, y_test = xdata, xdata, ydata, ydata
            fit_models, model_train_predictDF, train_stats = ModelTuningGrid(x_train=x_train, y_train=y_train, 
                ylabel=ylabel, sample_names=x_train.index.to_list(),  
                scaler=args.scaler, decom='none', method=fit_methods, scoring=args.scoring, random_state=args.random_state, n_jobs=args.threads)
        model_train_predictDF.to_csv(stats_path+'.train_predict.tsv', sep='\t', index=True, mode='w')
        train_stats.to_csv(stats_path+'.train_stats.tsv', sep='\t', index=True, mode='w')
        ROCforModels(fit_models=fit_models, x=x_train, y=y_train.to_list(), ylabel=ylabel, outpath=stats_path+'.train_data')
        # TODO: choose best model by accuracy/f1score/sensitivity/specificity/...
        
        if args.save is not None:
            model_path = args.save if args.save.endswith('.joblib') else args.save+'.joblib'
            os.makedirs(os.path.dirname(model_path), exist_ok=True)
            save_model = {'model':fit_models, 'ylabel':ylabel, 'features':features, 'refgenes':refgene_list, 'refgene_calculation':refgene_cal, 'transform':args.transform}
            joblib.dump(save_model, model_path)
    
    else: # only testing regardless of y_true_label
        load_model = joblib.load(args.model)
        fit_models = load_model['model']
        ylabel = load_model['ylabel']
        keep_features = load_model['features']
        count_matrix = MatrixTransform(count_matrix, method=load_model['transform'])
        if load_model['refgenes'] is not None:
            count_matrix = RelativeValuesByReferenceGenes(count_matrix, load_model['refgenes'], calculate_method=load_model['refgene_calculation'], log2level=True)
        
        x_test, y_test = PreprocessMatrixForTesting(matrix=count_matrix, ylabel=ylabel, group=sample_group, key='Group', features=keep_features)
        model_predictDF, test_stats = ModelTestingOnly(fit_models=fit_models, x_test=x_test, ylabel=ylabel, sample_names=x_test.index.to_list(), y_test=y_test.to_list())
        model_predictDF.to_csv(stats_path+'.test_predict.tsv', sep='\t', index=True, mode='w')
        test_stats.to_csv(stats_path+'.test_stats.tsv', sep='\t', index=True, mode='w')
        ROCforModels(fit_models=fit_models, x=x_test, y=y_test.to_list(), ylabel=ylabel, outpath=stats_path+'.test_data')
        
    return



if __name__ == '__main__':
    Main()

