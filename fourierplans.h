  //Setting up fftw_plans, in-place ffts:
  /*
  fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
			       fftw_complex *in, const int *inembed,
			       int istride, int idist,
			       fftw_complex *out, const int *onembed,
			       int ostride, int odist,
			       int sign, unsigned flags);
  */
  //  unsigned flags=FFTW_MEASURE;
  unsigned flags=( Vq > 1000 ? FFTW_PATIENT: FFTW_ESTIMATE);
  //unsigned flags=FFTW_PATIENT;

  logfile << "Making FFT plans" << endl;    
  FFTWCOMPLEX* A1_ptr  =reinterpret_cast<FFTWCOMPLEX*>(A1.start());
  FFTWCOMPLEX* A1r_ptr =reinterpret_cast<FFTWCOMPLEX*>(A1r.start());
  FFTWCOMPLEX* A2_ptr  =reinterpret_cast<FFTWCOMPLEX*>(A2.start());
  FFTWCOMPLEX* A2r_ptr =reinterpret_cast<FFTWCOMPLEX*>(A2r.start());
  FFTWCOMPLEX* B_ptr  =reinterpret_cast<FFTWCOMPLEX*>(B.start());
  FFTWCOMPLEX* Br_ptr =reinterpret_cast<FFTWCOMPLEX*>(Br.start());

  FFTWCOMPLEX* F1q_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F1q[0]);
  FFTWCOMPLEX* F1r_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F1r[0]);
  FFTWCOMPLEX* F2q_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F2q[0]);
  FFTWCOMPLEX* F2r_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F2r[0]);


#ifdef LONGDOUBLE
  A1q_to_A1r = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,-1,flags);
  A1r_to_A1q = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,+1,flags);
  A2q_to_A2r = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,-1,flags);
  A2r_to_A2q = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,+1,flags);
  Bq_to_Br   = fftwl_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,-1,flags);  
  Br_to_Bq   = fftwl_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,+1,flags);
  F1q_to_F1r = fftwl_plan_many_dft(dim,&dims[0],1,F1q_ptr,0,1,1,F1r_ptr,0,1,1,-1,flags);
  F1r_to_F1q = fftwl_plan_many_dft(dim,&dims[0],1,F1r_ptr,0,1,1,F1q_ptr,0,1,1,+1,flags);  
  F2q_to_F2r = fftwl_plan_many_dft(dim,&dims[0],1,F2q_ptr,0,1,1,F2r_ptr,0,1,1,-1,flags);
  F2r_to_F2q = fftwl_plan_many_dft(dim,&dims[0],1,F2r_ptr,0,1,1,F2q_ptr,0,1,1,+1,flags);  

#elif defined FLOAT
  A1q_to_A1r = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,-1,flags);
  A1r_to_A1q = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,+1,flags);
  A2q_to_A2r = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,-1,flags);
  A2r_to_A2q = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,+1,flags);
  Bq_to_Br   = fftwf_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,-1,flags);  
  Br_to_Bq   = fftwf_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,+1,flags);
  F1q_to_F1r = fftwf_plan_many_dft(dim,&dims[0],1,F1q_ptr,0,1,1,F1r_ptr,0,1,1,-1,flags);
  F1r_to_F1q = fftwf_plan_many_dft(dim,&dims[0],1,F1r_ptr,0,1,1,F1q_ptr,0,1,1,+1,flags);  
  F2q_to_F2r = fftwf_plan_many_dft(dim,&dims[0],1,F2q_ptr,0,1,1,F2r_ptr,0,1,1,-1,flags);
  F2r_to_F2q = fftwf_plan_many_dft(dim,&dims[0],1,F2r_ptr,0,1,1,F2q_ptr,0,1,1,+1,flags);  

#else
  A1q_to_A1r = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,-1,flags);
  A1r_to_A1q = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,+1,flags);
  A2q_to_A2r = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,-1,flags);
  A2r_to_A2q = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,+1,flags);
  Bq_to_Br   = fftw_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,-1,flags);  
  Br_to_Bq   = fftw_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,+1,flags);
  F1q_to_F1r = fftw_plan_many_dft(dim,&dims[0],1,F1q_ptr,0,1,0,F1r_ptr,0,1,0,-1,flags);
  F1r_to_F1q = fftw_plan_many_dft(dim,&dims[0],1,F1r_ptr,0,1,0,F1q_ptr,0,1,0,+1,flags);  
  F2q_to_F2r = fftw_plan_many_dft(dim,&dims[0],1,F2q_ptr,0,1,0,F2r_ptr,0,1,0,-1,flags);
  F2r_to_F2q = fftw_plan_many_dft(dim,&dims[0],1,F2r_ptr,0,1,0,F2q_ptr,0,1,0,+1,flags);  
#endif



  logfile << "Done making FFT plans" << endl;


