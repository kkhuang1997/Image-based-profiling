import feather
import pandas as pd
import sklearn.discriminant_analysis
import sklearn.ensemble
import sklearn.model_selection
import sklearn.svm
import sklearn.manifold
import sklearn.decomposition
import sklearn.calibration
import pickle
import numpy as np
import os


## load data
path = '../RDS data/features_BF_meta_drugs.feather'
df_feature_drug = feather.read_dataframe(path)
df = df_feature_drug
df_feature = df.iloc[:, 0:401].astype("float")
df_meta = df.iloc[:, 402:]

# lines = df_meta.Line.unique()
## drop "CK49P02" Line
lines = ['CK51P06', 'CK33P02', 'CK38P01', 'CK40P01', 'CK44P01']
for line in lines:
    acc_out_fn = os.path.join(
        "E:/classical feature extraction/table/drug_effect/SVM_Accuracies_Origin_{}.csv".format(line))
    profiles_mean_out_fn = os.path.join(
        "E:/classical feature extraction/table/drug_effect/SVM_Profiles_mean_Origin_{}.csv".format(line))
    profiles_std_out_fn = os.path.join(
        "E:/classical feature extraction/table/drug_effect/SVM_Profiles_std_Origin_{}.csv".format(line))

    features_line = df_feature.loc[df_meta.Line == line, :]
    metadata_line = df_meta.loc[df_meta.Line == line, :]

    features_dmso = features_line.loc[metadata_line.Drug == "DMSO", :]

    drugs = metadata_line.Drug.unique()
    drug_accuracies = pd.DataFrame(
        data=None,
        columns=["AUC_Mean", "AUC_Std", "Distance"],
        index=drugs)
    drug_profiles_mean = pd.DataFrame(
        data=None,
        columns=df_feature.columns.values,
        index=drugs)
    drug_profiles_std = pd.DataFrame(
        data=None,
        columns=df_feature.columns.values,
        index=drugs)

    for drug in drugs:
        print(line, "|", np.where(drugs == drug)[0][0] + 1, "of", len(drugs))
        features_drug = features_line.loc[metadata_line.Drug == drug, :]

        drug_dist = np.linalg.norm(  ## norm of Column Mean of treated and dmso
            features_drug.mean().values -
            features_dmso.mean().values)

        x = np.concatenate((features_dmso, features_drug), axis=0)
        y = np.repeat(("DMSO", "DRUG"), repeats=(
            len(features_dmso), len(features_drug)))

        # Train SVM
        x_param, x_train, y_param, y_train = \
            sklearn.model_selection.train_test_split(
                x, y, test_size=0.5)
        param_grid = [{"C": [5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2]}]
        svm_l1 = sklearn.model_selection.GridSearchCV(
            estimator=sklearn.linear_model.LogisticRegression(
                penalty="l2", dual=False),
            param_grid=param_grid, cv=5, scoring="roc_auc")
        svm_l1.fit(x_param, y_param == "DRUG")
        cv = list(sklearn.model_selection.KFold(
            n_splits=10, shuffle=True).split(X=x_train, y=y_train))
        clf = sklearn.calibration.CalibratedClassifierCV(
            base_estimator=sklearn.linear_model.LogisticRegression(
                penalty="l2", dual=False,
                C=svm_l1.best_estimator_.get_params()["C"]), cv=cv)
        clf.fit(X=x_train, y=y_train)

        aucs = []
        profiles = []
        for ii in range(len(cv)):
            try:
                train_i, test_i = cv[ii]
                estimator = clf.calibrated_classifiers_[ii]
                y_pred = estimator.predict_proba(X=x_train[test_i])
                aucs.append(sklearn.metrics.roc_auc_score(
                    y_true=(y_train[test_i] == "DRUG").astype(np.int),
                    y_score=y_pred[:, 1]))
                profiles.append(estimator.base_estimator.coef_)
            except ValueError:
                # continue
                pass

        drug_accuracies.loc[drug, :] = {
            "AUC_Mean": np.mean(aucs),
            "AUC_Std": np.std(aucs),
            "Distance": drug_dist}

        drug_profiles_mean.loc[drug, :] = np.mean(
            np.concatenate(profiles), axis=0)[None]

        drug_profiles_std.loc[drug, :] = np.std(
            np.concatenate(profiles), axis=0)[None]

    drug_accuracies.to_csv(acc_out_fn)
    drug_profiles_mean.to_csv(profiles_mean_out_fn)
    drug_profiles_std.to_csv(profiles_std_out_fn)
