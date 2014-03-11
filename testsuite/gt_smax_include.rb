Name "gt smax only index"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"   
  run_test "#{$bin}gt smax -esa Random.fna"
end

Name "gt smax absolute"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"   
  run_test "#{$bin}gt smax -esa Random.fna -absolute"
end

Name "gt smax diff absolute"
Keywords "gt_smax"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"
  run_test "#{$bin}gt smax -esa Random.fna -absolute"
  run_test "diff #{last_stdout} #{$testdata}Random.fna_smax_absolute.out"
end
