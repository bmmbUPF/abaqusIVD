C Initialization of variables
C ---------------------------

C Integer variables
      I = 0
      J = 0
      K = 0
      A = 0
      M = 0
      Ni = 0
      JJ = 0
      i_intPT = 0
      index = 0
      stat = 1

C--------------------------------------------------------------------------

C Initialize the shape functions and their local gradients
      N = ZERO
      dN = ZERO
      d2N = ZERO
      dNxi = ZERO
      d2Nxi = ZERO
      N0 = ZERO
      dN0 = ZERO
      d2N0 = ZERO
      dNC = ZERO
      d2NC = ZERO
      dN0C = ZERO
      d2N0C = ZERO
      d2Nmix = ZERO
      Nvec = ZERO
      DETMapJ = ZERO
      DETMapJ0 = ZERO
      DETMapJC = ZERO
      DETMapJ0C = ZERO
      Cip0 = ZERO
      Cip = ZERO
      whtG = ZERO
      WHT = ZERO

C Initialize the deformation gradient
      IdenM = ZERO
      IdenV = ZERO
      DFG = ZERO
      DFG_old = ZERO
      DFG_inv = ZERO
      DFG_invT = ZERO
      DETDFG = ZERO
      dDFGdx = ZERO

C--------------------------------------------------------------------------

C PROPS variables
      E = ZERO
      NU = ZERO
      D_WAT = ZERO
      ph_val_0 = ZERO
      CELL_rho_0 = ZERO
      CELL_viab_0 = ZERO
      ph_val_old = ZERO
      CELL_rho_old = ZERO
      CELL_viab_old = ZERO
      a_K = ZERO
      b_K = ZERO
      c_K = ZERO
      ALPHA_pH = ZERO
      GlucThres = ZERO
      pHThres = ZERO
      NF0h = ZERO

C I PROPS
      CNUMBER = 0
      ngSdv = 0
      nlSdv = 0
      NF0h_d = 0
      nD_d = 0
      PN_TIME = 0

C Transport variables
      NF = ZERO
      NF0 = ZERO
      pH_val = ZERO
      CELL_rho = ZERO
      CELL_viab = ZERO
      NDOFN = ZERO

C PN equation variables
      SP = ZERO

      TNF = ZERO
      TNF_PROT = ZERO

      IL1B = ZERO
      IL1B_PROT = ZERO

      AGG_IMN_ACC = ZERO
      AGG_IL1B _ACC= ZERO
      AGG_TNF_ACC = ZERO
      AGG_IL1B_TNF_ACC = ZERO
      
      COLI_IMN_ACC = ZERO
      COLI_IL1B_ACC = ZERO
      COLI_TNF_ACC = ZERO
      COLI_IL1B_TNF_ACC = ZERO

      COLII_IMN_ACC = ZERO
      COLII_IL1B_ACC = ZERO
      COLII_TNF_ACC = ZERO
      COLII_IL1B_TNF_ACC = ZERO

      MMP3_IMN_ACC = ZERO
      MMP3_IL1B_ACC = ZERO
      MMP3_TNF_ACC = ZERO
      MMP3_IL1B_TNF_ACC = ZERO

      ADM4_IMN_ACC = ZERO
      ADM4_IL1B_ACC = ZERO
      ADM4_TNF_ACC = ZERO
      ADM4_IL1B_TNF_ACC = ZERO

C--------------------------------------------------------------------------

C Initialize the DOFs
C Displacements
      U = ZERO
      dU = ZERO
      U_old = ZERO

C Solutes concentrations
C O2
      O2_CONC = ZERO
      dO2_CONC = ZERO
      O2_CONC_old = ZERO
      O2_CONC_0 = ZERO
C Lactate
      LACT_CONC = ZERO
      dLACT_CONC = ZERO
      LACT_CONC_old = ZERO
      LACT_CONC_0 = ZERO
C Glucose
      GLUC_CONC = ZERO
      dGLUC_CONC = ZERO
      GLUC_CONC_old = ZERO
      GLUC_CONC_0 = ZERO

C Current Coordinates
      COORDSC = ZERO
      Le = ZERO

C--------------------------------------------------------------------------

C O2 integration point variables
      O2_int = ZERO
      O2_int_old = ZERO
      dO2_dt = ZERO                              !temporal derivative of O2
      dO2dx = ZERO                               !spatial derivative of O2
      dO2dt = ZERO
      O2_int_0 = ZERO

C Lactate integration point variables
      LACT_int = ZERO
      LACT_int_old = ZERO
      dLACT_dt = ZERO                            !temporal derivative of lactate
      dLACTdx = ZERO                             !spatial derivative of lactate
      dLACTdt = ZERO
      LACT_int_0 = ZERO

C Glucose integration point variables
      GLUC_int = ZERO
      GLUC_int_old = ZERO
      dGLUC_dt = ZERO                            !temporal derivative of glucose
      dGLUCdx = ZERO                             !spatial derivative of glucose
      dGLUCdt = ZERO
      GLUC_int_0 = ZERO

C--------------------------------------------------------------------------

C Initialize the Diffusion and Reaction coefficients
      D_SOL = ZERO
      dD_O2dx = ZERO
      dD_LACTdx = ZERO
      dD_GLUCDx = ZERO
      dD_SOLdx_MAG = ZERO
      R_SOL = ZERO
      dR_SOL = ZERO

C--------------------------------------------------------------------------

C Initialize the residual and tangent matrices to zero.
      R_U = ZERO
      R_O2 = ZERO
      R_LACT = ZERO
      R_GLUC = ZERO
      K_U = ZERO
      K_O2 = ZERO
      K_O2LACT = ZERO
      K_O2GLUC = ZERO
      K_LACT = ZERO
      K_LACTO2 = ZERO
      K_LACTGLUC = ZERO
      K_GLUC = ZERO
      K_GLUCO2 = ZERO
      K_GLUCLACT = ZERO

C Initialize the residual and tangent factors to zero.
      RF_O2 = ZERO
      RF_LACT = ZERO
      RF_GLUC = ZERO
      TF_O2 = ZERO
      TF_LACT = ZERO
      TF_GLUC = ZERO