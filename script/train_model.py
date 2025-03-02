import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.ensemble import (BaggingRegressor, RandomForestRegressor,
                            ExtraTreesRegressor, GradientBoostingRegressor)
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR, LinearSVR
from sklearn.linear_model import Ridge
from xgboost import XGBRegressor
from script.load_descriptors import get_descriptors, f_des_std

def clean_data(X):
    """Clean data by replacing NaN and infinite values with 0."""
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return X

def train_and_evaluate():
    """Train and evaluate models using independent descriptors."""
    # Load data
    df = pd.read_csv('./dataset/all_input_data.csv')
    labels = np.array(df['Yield (%)'].to_list())
    labels_std = np.array(labels)/100

    # Load descriptors
    desc = get_descriptors()
    desc = clean_data(desc)  # Clean before standardization
    desc_std = f_des_std(desc)
    desc_std = clean_data(desc_std)  # Clean after standardization

    print(f"Descriptor shape: {desc_std.shape}")
    print(f"Number of samples: {len(labels_std)}")

    # Define models
    models = [
        BaggingRegressor(n_jobs=-1, random_state=42),
        DecisionTreeRegressor(random_state=42),
        ExtraTreesRegressor(n_jobs=-1, random_state=42),
        GradientBoostingRegressor(random_state=42),
        KNeighborsRegressor(n_jobs=-1),
        KernelRidge(),
        LinearSVR(max_iter=10000),
        RandomForestRegressor(n_jobs=-1, random_state=42),
        Ridge(random_state=42),
        SVR(),
        XGBRegressor(n_jobs=-1, random_state=42)
    ]
    model_names = ['BG', 'DT', 'ET', 'GB', 'KNR', 'KRR', 'LSVR', 'RF', 'Ridge', 'SVR', 'XGB']

    # Train and evaluate models
    kfold = KFold(n_splits=5, shuffle=True, random_state=42)
    performance_dict = {}

    for model, model_name in zip(models, model_names):
        print(f"\nTraining {model_name}...")
        all_r2 = []
        all_pearsr = []
        all_mae = []
        all_rmse = []
        repeat_pred = []
        repeat_test = []
        
        for train_index, val_index in kfold.split(desc_std):
            train_x, test_x = desc_std[train_index], desc_std[val_index]
            train_y, test_y = labels_std[train_index], labels_std[val_index]
            
            try:
                model.fit(train_x, train_y)
                test_pred = model.predict(test_x) * (labels.max() - labels.min()) + labels.min()
                test_y = test_y * (labels.max() - labels.min()) + labels.min()
                
                r2 = r2_score(test_y, test_pred)
                pearsr = pearsonr(test_y, test_pred)
                mae = mean_absolute_error(test_y, test_pred)
                rmse = (mean_squared_error(test_y, test_pred)) ** 0.5
                
                all_r2.append(r2)
                all_pearsr.append(pearsr[0])
                all_mae.append(mae)
                all_rmse.append(rmse)
                repeat_pred.append(test_pred)
                repeat_test.append(test_y)
            except Exception as e:
                print(f"Error training {model_name}: {str(e)}")
                continue
        
        if len(all_r2) > 0:
            performance_dict[model_name] = [
                np.mean(all_r2), np.mean(all_pearsr),
                np.mean(all_mae), np.mean(all_rmse),
                repeat_pred[np.argmax(all_r2)],
                repeat_test[np.argmax(all_r2)]
            ]
            
            print(f'Model: {model_name:>5}, R2: {np.mean(all_r2):.4f}, '
                  f'MAE: {np.mean(all_mae):.4f}, RMSE: {np.mean(all_rmse):.4f}, '
                  f'Pearson R: {np.mean(all_pearsr):.4f}')
        else:
            print(f"Model {model_name} failed to train on any fold")
    
    return performance_dict

if __name__ == '__main__':
    train_and_evaluate() 