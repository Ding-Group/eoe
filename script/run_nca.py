import pandas as pd
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.neighbors import KNeighborsClassifier

X_test = pd.read_csv('~/work/eoe/result/other/data_love.tsv', sep=' ', header=None)
X_test = X_test.to_numpy()

nca = NeighborhoodComponentsAnalysis(n_components=50, random_state=42)

for z in range(1, 11):
    z = str(z)

    name = '~/work/eoe/result/other/data_eoe_sample_' + z + '.tsv'
    X_train = pd.read_csv(name, sep=' ', header=None)
  
    name = '~/work/eoe/result/other/cell_type_eoe_num_sample_' + z + '.tsv'
    y_train = pd.read_csv(name, sep=' ', header=None)
  
    X_train = X_train.to_numpy()
    y_train = y_train.to_numpy()
  
    nca.fit(X_train, y_train)
  
    knn = KNeighborsClassifier(n_neighbors=11)
    knn.fit(nca.transform(X_train), y_train)
  
    print(knn.score(nca.transform(X_train), y_train))

    y_pred = knn.predict_proba(nca.transform(X_test))
  
    name = '~/work/eoe/result/other/cell_type_love_pred_prob_sample_' + z + '.tsv'
    pd.DataFrame(y_pred).to_csv(name,
                                header=False, index=False, sep='\t')
