# Python re-run of the LQR + integral-action design in matlab/main.m (same A, B, C1, Q_bar1, R).
# Usage: python lqr_tracking.py ../images/lqr_tracking.png
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.integrate import solve_ivp
l=0.235; Ix=Iy=8.11858e-5; Iz=6.12233e-5; m=0.506; g=9.8
A=np.zeros((12,12))
for i in range(0,12,2): A[i,i+1]=1
A[9,2]=g; A[11,0]=-g
B=np.zeros((12,4)); B[1,1]=l/Ix; B[3,2]=l/Iy; B[5,3]=l/Iz; B[7,0]=1/m
C1=np.zeros((4,12)); C1[0,4]=C1[1,6]=C1[2,8]=C1[3,10]=1   # psi, z, x, y
Ab=np.block([[A,np.zeros((12,4))],[-C1,np.zeros((4,4))]]); Bb=np.vstack([B,np.zeros((4,4))])
Q=10*np.eye(16); R=np.eye(4)
P=solve_continuous_are(Ab,Bb,Q,R); K=np.linalg.solve(R,Bb.T@P)
Acl=Ab-Bb@K
print('closed-loop stable:', np.all(np.linalg.eigvals(Acl).real<0), 'slowest pole:', np.max(np.linalg.eigvals(Acl).real).round(3))
def ref(t): return np.array([0.5*(t>=12), 1.0*(t>=1), 1.0*(t>=4), 1.0*(t>=8)])  # psi, z, x, y
def f(t,s): return Acl@s + np.concatenate([np.zeros(12), ref(t)])
T=np.linspace(0,25,2501); sol=solve_ivp(f,(0,25),np.zeros(16),t_eval=T,max_step=0.01)
Y=C1@sol.y[:12]; Rf=np.array([ref(t) for t in T]).T
labels=['yaw ψ (rad)','altitude z (m)','x (m)','y (m)']; cols=['#7c3aed','#0ea5e9','#f97316','#16a34a']
fig,ax=plt.subplots(figsize=(8,3.8))
for i in [1,2,3,0]:
    ax.plot(T,Rf[i],'--',color=cols[i],lw=1,alpha=.6); ax.plot(T,Y[i],color=cols[i],lw=2,label=labels[i])
ax.set_xlabel('time (s)'); ax.set_ylabel('output'); ax.set_title('LQR with integral action: step tracking (dashed = reference)')
ax.grid(alpha=.3); ax.legend(loc='lower right',fontsize=9); fig.tight_layout(); fig.savefig(sys.argv[1],dpi=130)
for i in range(4):
    yi,ri=Y[i],Rf[i]; tstep=T[np.argmax(ri>0)]; final=ri[-1]
    idx=np.where(np.abs(yi-final)>0.02*final)[0]; ts=T[idx[-1]]-tstep if len(idx) else 0
    print(f'{labels[i]}: overshoot {max(0,(yi.max()-final)/final*100):.1f}%  settling(2%) {ts:.2f}s  final err {abs(yi[-1]-final):.1e}')
