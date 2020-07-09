import pandas as pd
import numpy as np


# 超級平面状に一様に点を分布させる
# http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
# cell cycle scoreが小さくなるようなviewing axisを見つける。そのaxisは、phiとthetaで表す。
# ここで求めた、phi, thetaに従い、元PC達を全て、z軸周りに-phi, y軸周りに-theta回転させれば、
def uni_dist_on_unit_sphere(n=10000):
    goldenRatio = (1 + 5**0.5) / 2
    i = np.arange(0, n)
    phi = 2 * np.pi * i / goldenRatio
    theta = np.arccos(1 - 2*(i+0.5)/n)
    x, y, z = np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)
    return x, y, z, phi, theta # [3 * n]

# input: df [n_phases(=6) * n_pcs(=3)]
# そのPCに沿った各phase発現量の平均値
# 回転させるべき角度を返す
def find_optimal_viewing_axis(df_preprocessed_pca, df_phase, pc0, pc1, pc2, n=10000):
    min_cc_score = 10^4
    viewing_axis_idx = 0
    x, y, z, phi, theta = uni_dist_on_unit_sphere(n)
    
    df_preprocessed_pca_phase = pd.concat([df_preprocessed_pca, df_phase])
    df_preprocessed_pca_phase_mean = df_preprocessed_pca_phase.T.set_index("phase").astype(float).groupby("phase").mean()
    for idx, (i,j,k) in enumerate(zip(x,y,z)):
        # viewing axisの候補ベクトル
        candidate = np.array([i,j,k])
        ls = []
        for point in df_preprocessed_pca_phase_mean[[pc0, pc1, pc2]].iterrows():
            pnt =  np.array(point[1:][0])
            ls.append(np.dot(candidate,pnt))
        cc_score = np.var(ls)
        if cc_score < min_cc_score:
            min_cc_score = cc_score
            viewing_axis_idx = idx
    
    optimal_phi = phi[viewing_axis_idx] % (2*np.pi)
    optimal_theta = theta[viewing_axis_idx] % (np.pi)
    return optimal_phi, optimal_theta


# pc3のcell cycle scoreを小さくするように回転する
# pcの番号を渡す
# 回転させたいpcのdf
def rotation(df_preprocessed_pca, df_phase, pc0, pc1, pc2, n_components, n=10000):
    # 50PCたちを、z軸周りに-optimal_phiだけ回転して、y軸周りに-optimal_thetaだけ回転する
    # z軸周りに-optimal_phiだけ回転する用の行列
    optimal_phi, optimal_theta = find_optimal_viewing_axis(df_preprocessed_pca, df_phase, pc0, pc1, pc2, n)

    mtrx_phi = np.eye(n_components)
    mtrx_phi[pc0,pc0] = np.cos(optimal_phi)
    mtrx_phi[pc0,pc1] = np.sin(optimal_phi)
    mtrx_phi[pc1,pc0] = -np.sin(optimal_phi)
    mtrx_phi[pc1,pc1] = np.cos(optimal_phi)

    # y軸周りに-optimal_theta
    mtrx_theta = np.eye(50)
    mtrx_theta[pc0,pc0] = np.cos(optimal_theta)
    mtrx_theta[pc0,pc2] = -np.sin(optimal_theta)
    mtrx_theta[pc2,pc0] = np.sin(optimal_theta)
    mtrx_theta[pc2,pc2] = np.cos(optimal_theta)

    df_rotated = pd.DataFrame(mtrx_theta.dot(mtrx_phi.dot(df_preprocessed_pca)))
    df_rotated.columns = df_preprocessed_pca.columns
    return df_rotated

