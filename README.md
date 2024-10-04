# S-FAR
Program code of  regularization  and techniques: Including the S-FAR Methods, Debiasing Technique, and Thresholding Technique

"There are three folders : <**compare_with_CG**> folder, <**debiasing**> folder, and <**thresholding**> folder."

1.  <**compare_with_CG**> folder

>There is code for the S-FAR method with \theta =0.9, 1.5, 1.99,  Landweber method and  CG method. We provide a 1-dimensional example (see Example 1(2) in our paper) and a 2-dimensional example (see Example 3(2) in our paper).
>Please run <generate_noise.m> first, then <test_many.m>, and finally <plt.m>  for any example.

2.  <**debiasing**> folder
>There is code for the debiased S-FAR method with \theta =0.9, 1.5, 1.99, debiased  Landweber method. We provide a 1-dimensional example (see Example 1(2) in our paper) and a 2-dimensional example (see Example 3(2) in our paper).
>First, ensure that the data (data_EX12.mat or data_EX32.mat) from  <**compare_with_CG**> folder is successfully imported.
>Please run <test_many.m> first, then  <plt.m>  for any example.

3. <**thresholding**> folder

>There is code for the debiased S-FAR method with \theta =0.9, 1.5, 1.99, debiased  Landweber method. We provide a 1-dimensional example (see Example 1(4) in our paper) and a 2-dimensional example (see Example 3(3) in our paper).
>Please run <generate_noise.m> first, then <test_many.m>, and finally <plt.m>  for any example.





