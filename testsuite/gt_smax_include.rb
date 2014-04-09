Name "gt smax scan no index"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt smax -ii foo", :retval => 1 
end

Name "gt smax map no index"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt smax -ii foo -map", :retval => 1
end

Name "gt smax map diff absolute"
Keywords "gt_smax"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB -dna -suf -sds no" +
    " -lcp -des no -md5 no"
  run_test "#{$bin}gt smax -ii at1MB -absolute -l 50 -map"
  run_test "grep -v '^#' #{last_stdout}"
  run_test "diff #{last_stdout} #{$testdata}at1MB_smax_absolute_l50.out"
end

Name "gt smax map diff relative"
Keywords "gt_smax"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB -dna -suf -sds no" +
    " -lcp -des no -md5 no"
  run_test "#{$bin}gt smax -ii at1MB -l 50 -map"
  run_test "grep -v '^#' #{last_stdout}"
  run_test "diff #{last_stdout} #{$testdata}at1MB_smax_relative_l50.out"
end

Name "gt smax scan diff absolute"
Keywords "gt_smax"
 Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB -dna -suf -sds no" +
    " -lcp -des no -md5 no"
  run_test "#{$bin}gt smax -ii at1MB -absolute -l 50"
  run_test "grep -v '^#' #{last_stdout}"
  run_test "diff #{last_stdout} #{$testdata}at1MB_smax_absolute_l50.out"
end

Name "gt smax scan diff relative"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB -dna -suf -sds no" +
    " -lcp -des no -md5 no"
  run_test "#{$bin}gt smax -ii at1MB -l 50"
  run_test "grep -v '^#' #{last_stdout}"
  run_test "diff #{last_stdout} #{$testdata}at1MB_smax_relative_l50.out"
end
