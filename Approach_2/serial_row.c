#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>
#include<math.h>
#include<string.h>
#define constant 6.28318530718
#define CLK CLOCK_MONOTONIC

struct timespec diff(struct timespec start, struct timespec end){
	struct timespec temp;
	if((end.tv_nsec-start.tv_nsec)<0){
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}


typedef struct {
  unsigned char gs;
} PPMPixelGS;


typedef struct {
  int x, y;
  PPMPixelGS *data;
} PPMImageGS;

typedef struct{
	double real;
	double imag;
} Complex;

#define RGB_COMPONENT_COLOR 255
void writePPMGS(const char *filename, PPMImageGS *img);
static PPMImageGS *readPPMGS(const char *filename);

/*-----------------------------------convert to image complex arrays------------------------------*/
Complex** convert(PPMImageGS *im)
{
  int rows = im->x;
  int cols = im->y;
  int i,j,idx;
  Complex **arr = (Complex **)malloc(rows * sizeof(Complex *));
  for (i=0; i<rows; i++)
         arr[i] = (Complex *)malloc(cols * sizeof(Complex ));
  
  for(i=0;i<rows;i++)
  {
	for(j=0; j<cols; j++)
	{
		  idx = cols*i + j;
		  PPMPixelGS *temp = im->data + idx;
		  arr[i][j].real=(double)temp->gs;
		  arr[i][j].imag=0.0;
	}
  }
   return arr;
}
/*------------------------------look up table-----------------------------------------------*/
void twiddle( int n, double w[] )
{
  double arg;
  double aw;
  int i;
  int n2;
  const double pi = constant/2;

  n2 = n / 2;
  aw = 2.0 * pi / ( ( double ) n );

  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( double ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}

/*-------------------------------------FFT--------------------------------*/
void FFT(double *x, double *y, int n, double w[])
{

   int m,i,j,k,i2;
   double tx,ty;
   m=0;
   i= n;
   /* m = logN calculation*/
	while(i>0)
	{
		i/=2;
		m++;
	}
	m-=1;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n -1;i++)
   {
      if (i < j)
      {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }
   
 /*  FFT computation  */
   int mj, term_i, mi, j2, count2;
   double u1, u2, t1, t2;
   mj = 1;    //stride of j
   for(k=0;k<m;k++)
   {
	mi = 2*mj;	//stride of i
	term_i = n/mi;
		
	for(i=0; i<term_i; i++)
	{
		count2=0;
		for(j=i*mi;count2<mj;j++, count2++)
		{
			j%=(n-1);
			j2 = (j+ mj);

			int twiddle_index = count2*n/mi;
			u1 = w[twiddle_index*2+0];
			u2 = -w[twiddle_index*2+1];


			t1 = u1 * x[j2] - u2 * y[j2];
            		t2 = u1 * y[j2] + u2 * x[j2];
		        x[j2] = x[j] - t1;
		        y[j2] = y[j] - t2;
		        x[j] += t1;
		        y[j] += t2;
		}
		
	}
        mj = mj*2;
    }  
}

/*-----------------------------------2D FFT------------------------------*/
void FFT_2D(Complex **comp_in,int rows, int cols, double *w)
{
	
	
	int i,j;
	for(i=0;i<rows;i++)
	    {
		double x[rows];
		double y[rows];
		
	        for(j=0; j<cols; j++)
		{
			x[j]=comp_in[i][j].real;
			y[j]=comp_in[i][j].imag;
		}		
	        FFT(x,y,cols,w);
		
		for(j=0; j<cols; j++)
		{
			comp_in[i][j].real=x[j];
			comp_in[i][j].imag=y[j]	;
		}

	    }
}

/*-----------------------------------calculate the transpose------------------------------*/
Complex ** transpose(int N,Complex **comp_in){
	
  int blockrow, blockcolumn, i = 0, j = 0;
  int blocksize;
  blocksize = 16;
  
  Complex **arr = (Complex **)malloc(N * sizeof(Complex *));
  for (i=0; i<N; i++)
         arr[i] = (Complex *)malloc(N * sizeof(Complex ));
  //double start = omp_get_wtime();
  for (blockrow = 0; blockrow < N; blockrow += blocksize)
  {
	for (blockcolumn = 0; blockcolumn < N; blockcolumn += blocksize)
	{
		 for (i = blockrow; i < blockrow + blocksize; i++)
		 {
		          for (j = blockcolumn; j < blockcolumn + blocksize; j++)
			  {
				arr[i][j] = comp_in[j][i];
			  }
		 }
 	}
  }
  //double end = omp_get_wtime()-start;
  //printf("%.15f\n",end);
  return arr;
  
}
/*-----------------------------------convert complex arrays to image------------------------------*/

PPMImageGS * convert_comp_img(Complex **comp_in,int rows,int cols)
{
  int i,j;
  PPMImageGS *im2 = (PPMImageGS *) malloc(sizeof(PPMImageGS));
  im2->x = rows;
  im2->y = cols;
  im2->data = (PPMPixelGS *) malloc(rows*cols*sizeof(PPMPixelGS));
  double temp ;
  int idx;	
  	
  for(i=0;i<rows;i++)
    {
      for(j=0; j<cols; j++)
	{
	  idx = cols*i + j;
	  temp = sqrt(comp_in[i][j].real*comp_in[i][j].real + comp_in[i][j].imag*comp_in[i][j].imag);
	  PPMPixelGS *temp2 = im2->data + idx;
	  temp2->gs = floor(temp);
	}	
    }
   return im2;
}


/*-----------------------------------main function------------------------------*/
int main(int argc, char* argv[])
{

  struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
  clock_gettime(CLK, &start_e2e);

  int n = atoi(argv[1]);
  int p = atoi(argv[2]);
  //int run_id = atoi(argv[3]);
  char filename[30];
  strcpy(filename,"input/");
  strcat(filename,"gs_");
  strcat(filename,argv[1]);
  strcat(filename,".ppm");
  char *problem_name = "FFT";
  char *approach_name = "rows";

  PPMImageGS *image,*transformed_img;
  Complex **comp_in;
  image = readPPMGS(filename);
  int rows= image->x;
  int cols = image->y;
  double* w = ( double * ) malloc (rows* sizeof ( double ) );

  clock_gettime(CLK, &start_alg);

  comp_in = convert(image);

  twiddle(rows,w); 

  FFT_2D(comp_in, rows, cols,w);
	
  comp_in = transpose(rows,comp_in);
   
  FFT_2D(comp_in,rows,cols,w);
     
  comp_in = transpose(cols,comp_in);

  transformed_img = convert_comp_img(comp_in,rows,cols);

  clock_gettime(CLK, &end_alg);
  
  char out_file[30] ;
  strcpy(out_file,"output/");
  strcat(out_file,"gs_");
  strcat(out_file,argv[1]);
  strcat(out_file,"fft_serial.ppm");

  writePPMGS(out_file,transformed_img);
  free(w);

  clock_gettime(CLK, &end_e2e);

	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);
	
	//printf("%d,%d,%ld,%ld,%ld,%ld\n", n, p, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);
	
  
	printf("%s,%s,%d,%d,%ld,%ld,%ld,%ld\n",problem_name, approach_name, n, p, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);
	
  return 0;
}




/*-----------------------------------convert to image array to ppm------------------------------*/
void writePPMGS(const char *filename, PPMImageGS *img)
{
  FILE *fp;
  //open file for output
  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    exit(1);
  }

  //write the header file
  //image format
  fprintf(fp, "P5\n");

    

  //image size
  fprintf(fp, "%d %d\n",img->x,img->y);

  // rgb component depth
  fprintf(fp, "%d\n",RGB_COMPONENT_COLOR);

  // pixel data
  fwrite(img->data, img->x, img->y, fp);

  fclose(fp);
}


/*-----------------------------------convert image(ppm) to array------------------------------*/
static PPMImageGS *readPPMGS(const char *filename)
{
  char buff[16];
  PPMImageGS *img;
  FILE *fp;
  int c, rgb_comp_color;
  //open PPM file for reading
  fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    exit(1);
  }

  //read image format
  if (!fgets(buff, sizeof(buff), fp)) {
    perror(filename);
    exit(1);
  }

  //check the image format
  if (buff[0] != 'P' || buff[1] != '5') {
    fprintf(stderr, "Invalid image format (must be 'P5')\n");
    exit(1);
  }

  //alloc memory form image
  img = (PPMImageGS *)malloc(sizeof(PPMImageGS));
  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  //check for comments
  c = getc(fp);
  while (c == '#') {
    while (getc(fp) != '\n') ;
    c = getc(fp);
  }

  ungetc(c, fp);
  //read image size information
  if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
    fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
    exit(1);
  }

  //read rgb component
  if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
    fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
    exit(1);
  }

  //check rgb component depth
  if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
    fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
    exit(1);
  }

  while (fgetc(fp) != '\n') ;
  //memory allocation for pixel data
  img->data = (PPMPixelGS*)malloc(img->x * img->y * sizeof(PPMPixelGS));

  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  //read pixel data from file
  if (fread(img->data, img->x, img->y, fp) != img->y) {
    fprintf(stderr, "Error loading image '%s'\n", filename);
    exit(1);
  }


  fclose(fp);


  return img;
}
