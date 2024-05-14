import ee
import pandas as pd
import geopandas as gpd
import numpy as np
from sklearn.model_selection import train_test_split
import geemap
import smogn
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA

def construct_ee_model(df, label='spm', 
                       downsample=0,
                       n_estimators=100, test_size=0.30, q=5, 
                       log_transform=False, 
                       apply_smogn=False, smogn_settings={}):
    """ Construct GEE RF-regressor (ee.Classifier) trained and validated on provided DataFrame. """

    def _df_to_ee(df):
        """Convert a dataframe to a featurecollection with null geometries."""
        features = []
        for _, row in df.iterrows():
            properties = row.to_dict()
            feature = ee.Feature(None, properties)
            features.append(feature)
        return ee.FeatureCollection(features)

    def _ee_to_df(ee_fc, size=None, limit=50):
        """Convert ee.FeatureCollection to a dataframe using pagination."""
        def _ee_subset_to_df(ee_fc, offset):
            ee_fc_temp = ee.FeatureCollection(ee_fc.toList(limit, offset))
            df_temp = geemap.ee_to_df(ee_fc_temp)
            return df_temp
        if size == None: size = ee_fc.size().getInfo()
        offsets = range(0, size, limit)
        dfs = [_ee_subset_to_df(ee_fc, offset) for offset in offsets]
        df = pd.concat(dfs, ignore_index=True)
        return df

    # setup train/test sets
    features = df.copy()
    
    if downsample:
        features_sz = features.shape[0]
        features = features.sample(round(features.shape[0]*downsample))
        print(f'Downsampling by {downsample*100:0.2f}% from {features_sz} to {features.shape[0]} samples.')

    labels = np.array(features[label])
    features['bin'] = pd.qcut(features[label], q=q) # quantile binning
    #features['bin'] = pd.cut(features[label], bins=q) # intervals binning
    bins = np.array(features['bin'])

    df_train, df_test, train_labels, test_labels = train_test_split(
        features, labels, 
        test_size=test_size, 
        stratify=bins, 
        random_state=42)
    df_train = df_train.drop(columns=['bin']).reset_index().drop(columns=['index'])
    df_test = df_test.drop(columns=['bin']).reset_index().drop(columns=['index'])

    if apply_smogn:
        print('Applying SMOGN to create balanced training set.')
        train_sz_bef = df_train.shape[0]
        df_train_smogn = smogn.smoter(
            data=df_train,
            y=label,
            **smogn_settings
        )
        train_sz_after = df_train.shape[0]
        traing_sz_d = (train_sz_after-train_sz_bef)/train_sz_bef*100
        print(f'Training set size from {train_sz_bef} to {train_sz_after} ({traing_sz_d:0.1f}%).')
        # df_train[label].plot.kde(label='Original')
        # df_train_smogn[label].plot.kde(label=f'SMOGN ({smogn})')
        df_train = df_train_smogn
        
    if log_transform: # log transform
        df_train['spm'] = np.log1p(df_train['spm'])
        df_test['spm'] = np.log1p(df_test['spm'])

    # create and train rf-regressor
    ee_fc_train = _df_to_ee(df_train)
    ee_fc_test = _df_to_ee(df_test)
    ee_rf_model = ee.Classifier.smileRandomForest(**{
        'numberOfTrees': n_estimators,
        'seed': 42
        }) \
        .setOutputMode('REGRESSION') \
        .train(**{
            'features': ee_fc_train,
            'classProperty': label,
            'inputProperties': list(df_train.columns)
        })
    ee_fc_predicted = ee_fc_test.classify(ee_rf_model, outputName=f'{label}_predicted')

    df_predicted = _ee_to_df(ee_fc_predicted, size=df_test.shape[0])

    if log_transform: # inverse log transform
        df_train[label] = np.expm1(df_train[label])
        df_test[label] = np.expm1(df_test[label])
        df_predicted[label] = np.expm1(df_predicted[label])
        df_predicted[f'{label}_predicted'] = np.expm1(df_predicted[f'{label}_predicted'])

    return(df_train, df_test, df_predicted, ee_rf_model)

def construct_local_model(df, label='spm', 
                       n_estimators=100, n_jobs=-1,
                       test_size=0.30, q=5, 
                       log_transform=False, 
                       apply_pca=False, n_components=2,
                       apply_smogn=False, smogn_settings={}):
    """ Construct a local RF-regressor (sklearn.ensemble.RandomForestRegressor) 
        trained and validated on provided DataFrame. """

    # setup train/test sets
    features = df.copy()
    labels = np.array(features[label])

    if apply_pca:
        pca = PCA(n_components=n_components)
        features = features.drop(columns=[label])
        pca.fit(features)
        features_pca = pd.DataFrame(pca.transform(features))
        features_pca = features_pca.rename(columns=lambda x: f'comp_{x}')
        features_pca[label] = labels
        features = features_pca
    
    # stratification bins
    features['bin'] = pd.qcut(features[label], q=q) # quantile binning
    #features['bin'] = pd.cut(features[label], bins=q) # intervals binning
    bins = np.array(features['bin'])

    df_train, df_test, train_labels, test_labels = train_test_split(
        features, labels, 
        test_size=test_size, 
        stratify=bins, 
        random_state=42)
    df_train = df_train.drop(columns=['bin']).reset_index().drop(columns=['index'])
    df_test = df_test.drop(columns=['bin']).reset_index().drop(columns=['index'])

    if apply_smogn:
        print('Applying SMOGN to create balanced training set.')
        train_sz_bef = df_train.shape[0]
        df_train_smogn = smogn.smoter(
            data=df_train,
            y=label,
            **smogn_settings
        )
        train_sz_after = df_train.shape[0]
        traing_sz_d = (train_sz_after-train_sz_bef)/train_sz_bef*100
        print(f'Training set size from {train_sz_bef} to {train_sz_after} ({traing_sz_d:0.1f}%).')
        # df_train[label].plot.kde(label='Original')
        # df_train_smogn[label].plot.kde(label=f'SMOGN ({smogn})')
        df_train = df_train_smogn
        train_labels = df_train.spm
        
    if log_transform: # log transform
        df_train['spm'] = np.log1p(df_train['spm'])
        df_test['spm'] = np.log1p(df_test['spm'])

    rf = RandomForestRegressor(n_estimators=n_estimators, random_state=42, n_jobs=n_jobs, verbose=True)
    rf.fit(df_train.drop(columns=[label]), train_labels)
    df_predicted = df_test.copy()
    df_predicted[f'{label}_predicted'] = rf.predict(df_test.drop(columns=[label]))

    if log_transform: # inverse log transform
        df_train[label] = np.expm1(df_train[label])
        df_test[label] = np.expm1(df_test[label])
        df_predicted[label] = np.expm1(df_predicted[label])
        df_predicted[f'{label}_predicted'] = np.expm1(df_predicted[f'{label}_predicted'])

    return(df_train, df_test, df_predicted, rf)