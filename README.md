# Stability-OPF

This code supplements the material presented in the following publications:

1. M. Bazrafshan, N. Gatsis, A. F. Taha, and J. A. Taylor, “Augmenting the optimal power flow for stability,” in Proc. 55th IEEE Conference on Decision and Control (CDC), Las Vegas, NV, 2016, pp. 4104-4109.

2. M. Bazrafshan, N. Gatsis, A. Taha,  and J. A. Taylor, “Coupling load-following stability with OPF,” IEEE Trans. Smart Grid, 2017, [Submitted]. 

# Usage
1. Open main_script.m. 
2. You can choose either of the following casefiles. The casefiles are located in Stability-OPF/casefiles
  * case9wmac_con
  * case14wmac_con
  * case39wmac_con
3. Select `stcontrol='OPF'` or `stcontrol='LQR-OPF'`. 
4. Select `lfcontrol='LQR'`. 
5. Choose the coupling parameter `alpha=0` (or any other value between 0 and 1). 
6. Run the command `workflow(casefile, stcontrol, lfcontrol,alpha)`. 


