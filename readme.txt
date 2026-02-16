プログラムの内容

bridge.msh
・Gmshで作成した橋形の構造物にランダムメッシュを生成させたもの

MSH_to_csv.c
・bridge.mshの不要なノードとエッジの情報を取り除いてcsvファイルにするコード

0122.c
・応力依存成長則の３次元版を使って変位させるコード

mesh_3_deformed_step_~~.csv
・計算ステップ100回ごとにその時点での変位の状況を記録したcsvファイル

edge_radius_history.csv(.png)
・計算ステップごとの全てのエッジの太さを示したグラフとcsvファイル

New_mesh_ope.c
・3Dプリンタようにエッジに補正を加えた結果とエッジを均一化させた結果のcsvファイルを出力するコード
	出力ファイル１：mesh_raw_2.csv(エッジ補正後の応力依存の構造物)
	出力ファイル２：mesh_uni_2.csv(エッジを均一化させた構造物)

CtoSCAD.c
・OPENSCADにそのままコピペできるtxtファイルを出力するコード

