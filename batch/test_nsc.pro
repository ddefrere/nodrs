PRO TEST_NSC
  n = 10
  date = '170406'
  FOR i=0, n-1 DO REDUCE_HOSTS2, date, /reduce, /skip_red, /skip_flx, /renew, /fit, /calib_null
END