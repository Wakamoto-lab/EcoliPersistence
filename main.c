//
//  main.c
//  EcoliPersistence
//
//  Created on 2020/07/07
//  Last modified on 2022/05/03
//  Contributors: Miki Umetani, Yuichi Wakamoto
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>

#define GNUPLOT_PATH "/opt/homebrew/bin/gnuplot"
//#define GNUPLOT_PATH "gnuplot -persist"

int FlNormIndex = 0;
int TLInterval = 3;
//double camera_scale = 0.064428; //For ORCA-R2, x100 (NA 1.45) obj.
double camera_scale = 0.086806; //For ORCA-Flash, x100 (NA 1.45) obj.

typedef struct
{
  int index;
  double area;
  double fl;
  double xm;
  double ym;
  int slice;
  int lastIndex;
  int prevCell;
  int isConnected;
  int isNode;
  int nextCell[2];
  int nextCellNum;
  double backLevel;
  int deadROI;
} ROI;

int FileExistCheck(char fname[]);
int maxSliceReturn(ROI roi[], int roiCount);
int MaxSliceReturnFromAllFile(void);

long FileLineNumberReturn(char fname[]);
long RoiCounterInFile(char fname[]);

double Correl(double mean_x, double mean_y, double SD_x, double SD_y, double xy_mean);
double CorrelFromFile(char fname[]);
double MinValueReturnFile(char fname[]);
double ReturnLastTenFlMean(ROI *roi, int j, long maxSlice);
double ReturnLastTenPointElRate(ROI *roi, int j, long maxSlice);

void AddNumberExtToFileName(char fname[], int rank, char extention[], long number);
void ChangeExtentionTo(char fname[], char ext[]);
void ColumnInitialize_double(double column[], long clnum);
void ColumnInitialize_int(int column[], long clnum);
void DataExtract(char fname[], ROI *roi, int *roiCount);
void ElRateRpoSCorrelation(void);
void eps_output(FILE *gp, char fname[]);
void gpsave(FILE *gp, char fname[]);
void gpStandardFormat(FILE *gp);
void InputNonNegativeIntValue(int *value, char sentence[]);
void InputPositiveDoubleValue(double *value, char sentence[]);
void LastFlDistribution(void);
void LastTenPointElongationRateDist(void);
void MeanVarAdd(double *mean, double *var, double value);
void MeanVarCalcFile(char fname[], double *mean, double *var);
void MeanVarFinal(double *mean, double *var, double counter);
void NormalizingFlValue(ROI roi[], int roiCount);
void NumberOrder(long number, int *order);
void NumberStringReturn(char numberstring[], int rank, long number);
void PersistenceAnalysis(void);
void plot_histogram_width_specified(char filename[], char xlabel[], char ylabel[], float fill_level, int eps_output_index, int datanum, double binSize, double binRatio);
void plot_linespoints(char filename[], char plottype[], char xlabel[], char ylabel[], int xlogscale, int ylogscale, int eps_output_index, int datanum);
void RoiNextCellAssign(ROI roi[], int i);

int main(){
  int index;

  if(camera_scale < 0.0645 && camera_scale > 0.0644)
    printf("Notice: Scale condition is adjusted to 'Hamamatsu ORCA-R2'\n");
  if(camera_scale < 0.0869 && camera_scale > 0.0868)
    printf("Notice: Scale condition is adjusted to 'Hamamatsu ORCA-Flash'\n");

  PersistenceAnalysis();

  return 0;
}

int FileExistCheck(char fname[])
{
  int returnValue;
  FILE *fin;
  fin = fopen(fname, "r");
  if(fin == NULL){
    returnValue = 0;
  }else{
    returnValue = 1;
  }
  fclose(fin);
  return returnValue;
}

int maxSliceReturn(ROI roi[], int roiCount)
{
  int maxSlice = 0;
  int i;
    
  for(i=0; i<roiCount; i++){
    if(roi[i].isConnected == 1){
      if(roi[i].slice > maxSlice){
        maxSlice = roi[i].slice;
      }
    }
  }
  return maxSlice;
}

int MaxSliceReturnFromAllFile(void)
{
  int i;
  //ROI roi[ROIMAX];
  ROI *roi;
  int roiCount;
  char fname[1000];
  int maxSliceThis;
  int maxSlice = 0;
  FILE *fin;

  for(i=0; i<10000; i++){
    strcpy(fname, "Data/Results");
    AddNumberExtToFileName(fname, 4, ".xls", i);
    fin = fopen(fname, "r");
    if(fin != NULL){
      fclose(fin);
      roi = (ROI *)malloc(sizeof(ROI)*RoiCounterInFile(fname));
      if(roi == NULL){
    printf("roi pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
      }
      DataExtract(fname, roi, &roiCount);
      maxSliceThis = maxSliceReturn(roi, roiCount);
      if(maxSliceThis > maxSlice) maxSlice = maxSliceThis;
      free(roi);
    }else{
      fclose(fin);
    }
  }
  return maxSlice;
}

long FileLineNumberReturn(char fname[])
{
  long counter = 0;
  FILE *fin;
  char buff[256];

  fin = fopen(fname, "r");
  if(fin == NULL){
    printf("ERROR: '%s' not exist!\n", fname);
  }else{
    while(fgets(buff, 256, fin)!=NULL){
      counter++;
    }
  }
  fclose(fin);

  return counter;
}

long RoiCounterInFile(char fname[])
{
  long i;
  char temp_label[1000];
  double background;
  ROI roiTemp;
  FILE *fin = fopen(fname, "r");

  for(i=0; i<9; i++){
    fscanf(fin, "%s", temp_label);
  }

  i=0;
  while(fscanf(fin, "%d", &roiTemp.index)!=EOF){
    fscanf(fin, "%s\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf", temp_label, &roiTemp.area, &roiTemp.fl, &roiTemp.xm, &roiTemp.ym, &roiTemp.slice, &roiTemp.lastIndex, &roiTemp.prevCell, &background);
    i++;
  }

  fclose(fin);

  return i;
}

double Correl(double mean_x, double mean_y, double SD_x, double SD_y, double xy_mean)
{
  return (xy_mean-mean_x*mean_y)/SD_x/SD_y;
}

double CorrelFromFile(char fname[])
{
  double counter=0.0;
  double x, y;
  double mean_x=0.0, mean_y=0.0, SD_x=0.0, SD_y=0.0, xy_mean=0.0;
  FILE *fin;

  fin = fopen(fname, "r");

  while(fscanf(fin, "%lf", &x)!=EOF){
    fscanf(fin, "%lf", &y);
    MeanVarAdd(&mean_x, &SD_x, x);
    MeanVarAdd(&mean_y, &SD_y, y);
    xy_mean = xy_mean + x*y;
    counter = counter + 1.0;
  }

  MeanVarFinal(&mean_x, &SD_x, counter);
  MeanVarFinal(&mean_y, &SD_y, counter);
  xy_mean = xy_mean/counter;

  fclose(fin);

  return Correl(mean_x, mean_y, pow(SD_x, 0.5), pow(SD_y, 0.5), xy_mean);
}

double MinValueReturnFile(char fname[])
{
  double valueTemp;
  double valueMin;
  FILE *fin;

  fin = fopen(fname, "r");

  if(fscanf(fin, "%lf", &valueMin)!=EOF){ //assigning initial value to min.
    while(fscanf(fin, "%lf", &valueTemp)!=EOF){
      if(valueTemp < valueMin){
    valueMin = valueTemp;
      }
    }
  }

  fclose(fin);
  return valueMin;
}

double ReturnLastTenFlMean(ROI *roi, int j, long maxSlice)
{
  double sum = 0.0, counter = 0.0;
  int i;

  i=j;
  while(roi[i].slice >= maxSlice - 10){
    sum = sum + roi[i].fl;
    counter = counter + 1.0;
    i = roi[i].prevCell;
  }

  return sum/counter;
}

double ReturnLastTenPointElRate(ROI *roi, int j, long maxSlice)
{
  int i;
  double sum = 0.0;
  double counter = 0.0;

  i = roi[j].prevCell;
  while(roi[i].slice >= maxSlice - 10){
    if(roi[i].isNode == 0){
      sum = sum
    + 1.0/TLInterval/(roi[j].slice-roi[i].slice)*log(roi[j].area/roi[i].area);
      counter = counter + 1.0;
    }
    j = i;
    i = roi[i].prevCell;
  }
  if(counter < 0.5){
    printf("Some lineage is strange!\n");
    exit(1);
  }

  return sum/counter;
}

void AddNumberExtToFileName(char fname[], int rank, char extention[], long number)
{
  int order=1;
  char fname_temp[100];
  char numberstring[100];

  NumberOrder(number, &order);//Determination of order
  if(rank<order){
    printf("Number is too large in 'AddNumberExtToFileName'\n");
    exit(EXIT_FAILURE);
  }

  NumberStringReturn(numberstring, rank, number);
  sprintf(fname_temp, "%s%s%s", fname, numberstring, extention);
  strcpy(fname, fname_temp);
}

void ChangeExtentionTo(char fname[], char ext[])
{
  int i;
  char fname_temp[1000];

  for(i=0; i<strlen(fname); i++){
    if(fname[i] == '.'){
      fname[i] = '\0';
    }
  }

  sprintf(fname_temp, "%s%s", fname, ext);
  strcpy(fname, fname_temp);

  //printf("File name changed to '%s'\n", fname);
}

void ColumnInitialize_double(double column[], long clnum)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=0.0;
  }
}

void ColumnInitialize_int(int column[], long clnum)
{
  long i;
  for(i=0;i<clnum;i++){
    column[i]=0;
  }
}

void DataExtract(char fname[], ROI *roi, int *roiCount)
{
  int i, j;
  char temp_label[1000];
  double background;
  ROI roiTemp;
  int *nodeCounter;
  int *usedIndex;
  FILE *fin;

  fin = fopen(fname, "r");

  for(i=0; i<9; i++){
    fscanf(fin, "%s", temp_label);
  }

  i=0;
  while(fscanf(fin, "%d", &roiTemp.index)!=EOF){
    fscanf(fin, "%s\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf",
       temp_label, &roiTemp.area, &roiTemp.fl, &roiTemp.xm, &roiTemp.ym,
       &roiTemp.slice, &roiTemp.lastIndex, &roiTemp.prevCell, &background);
    roiTemp.index = roiTemp.index - 1;
    roiTemp.prevCell = roiTemp.prevCell - 1;
    roiTemp.fl = roiTemp.fl - background;

    roiTemp.area = roiTemp.area * camera_scale * camera_scale;
    roi[i] = roiTemp;

    i++;
  }

  (*roiCount) = i;

  //Assign 'isConnected' and 'nextCellNum' values
  for(i=0; i<(*roiCount); i++){
    roi[i].isConnected = 0; //Assign 0 initially
    roi[i].nextCellNum = 0;
  }

  for(i=0; i<(*roiCount); i++){
    if(roi[i].lastIndex == 1){
      roi[i].isConnected = 1;
      j = roi[i].prevCell;
      while(1){
        if(j < 0){
         break;
        }
        roi[j].isConnected = 1;
        j = roi[j].prevCell;
      }
    }
  }

  usedIndex = (int *)malloc(sizeof(int)*(*roiCount));
  if(usedIndex == NULL){
    printf("usedIndex pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
  }
  ColumnInitialize_int(usedIndex, (*roiCount));
  for(i=0; i<(*roiCount); i++){
    if(roi[i].lastIndex == 1){
      RoiNextCellAssign(roi, i);
      usedIndex[i] = 1;
      j = roi[i].prevCell;
      while(1){
        if(usedIndex[j] == 1 || j < 0){
          break;
        }
        RoiNextCellAssign(roi, j);
        usedIndex[j] = 1;
        j = roi[j].prevCell;
      }
    }
  }
  free(usedIndex);

  //Assign 'isNode' values
  //printf("isNode value is being checked\n");
  nodeCounter = (int *)(malloc(sizeof(int)*(*roiCount)));
  if(nodeCounter == NULL){
    printf("'nodeCounter' couldn't be assigned in 'DataExtract()'!\n");
    exit(EXIT_FAILURE);
  }
  ColumnInitialize_int(nodeCounter, (*roiCount));
  for(i=0; i<(*roiCount); i++){
    if(roi[i].isConnected == 1 && roi[i].prevCell >= 0){
      nodeCounter[roi[i].prevCell]++;
    }
  }
  for(i=0; i<(*roiCount); i++){
    if(nodeCounter[i] == 2){
      roi[i].isNode = 1;
    }else if(nodeCounter[i] == 1){
      roi[i].isNode = 0;
    }else if(nodeCounter[i] > 2){
      printf("More than 2 lineages are converging to roi[%d]! when assigning node\n", i);
      exit(EXIT_FAILURE);
    }else{
      roi[i].isNode = 0;
    }
  }

  if(FlNormIndex == 1 && FileExistCheck("MeanFlForNormalization.dat") == 1){
    NormalizingFlValue(roi, (*roiCount));
  }

  //Assign roi[].deadROI = 0 initially, this is used in the chamber
  //array analysis
  for(i=0; i<(*roiCount); i++){
    roi[i].deadROI = 0;
  }

  free(nodeCounter);
  fclose(fin);
}

void ElRateRpoSCorrelation(void)
{
  long maxSliceAll;
  ROI *roi;
  int roiCount;
  long i,j;
  char fname[1000];
  char fname_out[] = "ElRateRpoSCorr.dat";
  //char fname_hist[] = "LastFlHist.dat";
  double corr;
  FILE *fin, *fout;
  
  fout = fopen(fname_out, "w");

  maxSliceAll = MaxSliceReturnFromAllFile();

  for(i=0; i<10000; i++){
    strcpy(fname, "Data/Results");
    AddNumberExtToFileName(fname, 4, ".xls", i);
    fin = fopen(fname, "r");
    if(fin != NULL){
      fclose(fin);

      roi = (ROI *)malloc(sizeof(ROI)*RoiCounterInFile(fname));
      if(roi == NULL){
    printf("roi pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
      }

      DataExtract(fname, roi, &roiCount);
      for(j=0; j<roiCount; j++){
    if(roi[j].slice == maxSliceAll && roi[j].lastIndex == 1){
      fprintf(fout, "%lf\t%lf\n",
          ReturnLastTenPointElRate(roi, j, maxSliceAll),
          ReturnLastTenFlMean(roi, j, maxSliceAll) );
    }
      }
      free(roi);
    }else{
      fclose(fin);
    }
  }
  fclose(fout);

  corr = CorrelFromFile(fname_out);
  printf("Correlation coefficient: %lf +/- %lf\n",
     corr, pow((1.0-corr*corr)/(FileLineNumberReturn(fname_out)-2), 0.5));

  plot_linespoints(fname_out, "p", "Elongation rate /min-1", "RpoS-mCherry fluorescence intensity /a.u.", 0, 0, 1, 1);

}

void eps_output(FILE *gp, char fname[])
{
  char fname_eps[256];

  strcpy(fname_eps, fname);
  ChangeExtentionTo(fname_eps, ".eps");

  fprintf(gp, "set term postscript eps enhanced color\n");
  fprintf(gp, "set output '%s'\n", fname_eps);
  fprintf(gp, "replot\n");
  //fprintf(gp, "set term x11\n");
  fprintf(gp, "set term qt\n");
  fprintf(gp, "replot\n");
}

void gpsave(FILE *gp, char fname[])
{
  char fname_eps[256];

  strcpy(fname_eps, fname);
  ChangeExtentionTo(fname_eps, ".plt");

  fprintf(gp, "save '%s'\n", fname_eps);
}

void gpStandardFormat(FILE *gp)
{
  fprintf(gp, "set size ratio 0.618\n");
  fprintf(gp, "set xtics out\n");
  fprintf(gp, "set xtics nomirror\n");
  fprintf(gp, "set ytics out\n");
  fprintf(gp, "set ytics nomirror\n");
  fprintf(gp, "set border 3\n");
}

void InputNonNegativeIntValue(int *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%d", value);
    if((*value)<0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}

void InputPositiveDoubleValue(double *value, char sentence[])
{
  while(1){
    printf("%s: ", sentence);
    scanf("%lf", value);
    if((*value)<=0.0){
      printf("Enter positive value!\n");
    }else{
      break;
    }
  }
}

void LastFlDistribution(void)
{
  long maxSliceAll;
  ROI *roi;
  int roiCount;
  long i,j;
  char fname[1000];
  double mean, var;
  double min, shift, binWidth;
  char fname_out[] = "LastFlData.dat";
  char fname_hist[] = "LastFlHist.dat";
  double column[1000];
  int pos;
  long dataNum;
  double value;
  FILE *fin, *fout;
  
  fout = fopen(fname_out, "w");

  maxSliceAll = MaxSliceReturnFromAllFile();

  for(i=0; i<10000; i++){
    strcpy(fname, "Data/Results");
    AddNumberExtToFileName(fname, 4, ".xls", i);
    fin = fopen(fname, "r");
    if(fin != NULL){
      fclose(fin);

      roi = (ROI *)malloc(sizeof(ROI)*RoiCounterInFile(fname));
      if(roi == NULL){
    printf("roi pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
      }

      DataExtract(fname, roi, &roiCount);
      for(j=0; j<roiCount; j++){
    if(roi[j].slice == maxSliceAll && roi[j].lastIndex == 1){
      fprintf(fout, "%lf\n", ReturnLastTenFlMean(roi, j, maxSliceAll));
    }
      }
      free(roi);
    }else{
      fclose(fin);
    }
  }
  fclose(fout);

  MeanVarCalcFile(fname_out, &mean, &var);
  printf("Mean: %lf min-1\n", mean);
  printf("SD: %lf min-1\n", pow(var, 0.5));
  printf("CV: %lf\n", pow(var, 0.5)/mean);
  dataNum = FileLineNumberReturn(fname_out);
  printf("N: %ld\n", dataNum);
  min = MinValueReturnFile(fname_out);
  printf("Min value: %lf\n", min);

  InputPositiveDoubleValue(&binWidth, "Enter bin width (75 in the manuscript)");

  //To make the histogram even for the negative values
  if(min < 0.0){
    shift = 0.0;
    while(shift > min){
      shift = shift - binWidth;
    }
  }else{
    shift = 0.0;
  }

  ColumnInitialize_double(column, 1000);
  fin = fopen(fname_out, "r");
  while(fscanf(fin, "%lf", &value)!=EOF){
    pos = (int)((value - shift)/binWidth);
    column[pos] = column[pos] + 1.0;
  }
  fclose(fin);

  for(i=0; i<1000; i++){
    column[i] = column[i]/dataNum;
  }

  fout = fopen(fname_hist, "w");
  for(i=0; i<1000; i++){
    if(column[i] > 0.0){
      fprintf(fout, "%lf\t%lf\n", shift + (i+0.5)*binWidth, column[i]);
    }
  }
  fclose(fout);

  plot_histogram_width_specified(fname_hist, "Fluorescence intensity /a.u.", "Frequency",
                 0.5, 1, 1, binWidth, 1.0);

}


void LastTenPointElongationRateDist(void)
{
  long maxSliceAll;
  ROI *roi;
  int roiCount;
  long i,j;
  char fname[1000];
  double mean, var;
  double min, shift, binWidth;
  char fname_out[] = "LastTenElRate.dat";
  char fname_hist[] = "LastTenElRateHist.dat";
  double column[1000];
  int pos;
  long dataNum;
  double value;
  FILE *fin, *fout;
  
  fout = fopen(fname_out, "w");

  maxSliceAll = MaxSliceReturnFromAllFile();

  for(i=0; i<10000; i++){
    strcpy(fname, "Data/Results");
    AddNumberExtToFileName(fname, 4, ".xls", i);
    fin = fopen(fname, "r");
    if(fin != NULL){
      fclose(fin);

      roi = (ROI *)malloc(sizeof(ROI)*RoiCounterInFile(fname));
      if(roi == NULL){
    printf("roi pointer couldn't be assigned!\n");
    exit(EXIT_FAILURE);
      }

      DataExtract(fname, roi, &roiCount);
      for(j=0; j<roiCount; j++){
    if(roi[j].slice == maxSliceAll && roi[j].lastIndex == 1){
      fprintf(fout, "%lf\n", ReturnLastTenPointElRate(roi, j, maxSliceAll));
    }
      }
      free(roi);
    }else{
      fclose(fin);
    }
  }
  fclose(fout);

  MeanVarCalcFile(fname_out, &mean, &var);
  printf("Mean: %lf min-1\n", mean);
  printf("SD: %lf min-1\n", pow(var, 0.5));
  printf("CV: %lf\n", pow(var, 0.5)/mean);
  dataNum = FileLineNumberReturn(fname_out);
  printf("N: %ld\n", dataNum);
  min = MinValueReturnFile(fname_out);
  printf("Min value: %lf\n", min);

  InputPositiveDoubleValue(&binWidth, "Enter bin width (0.002 in the manuscript)");

  //To make the histogram even for the negative values
  if(min < 0.0){
    shift = 0.0;
    while(shift > min){
      shift = shift - binWidth;
    }
  }else{
    shift = 0.0;
  }

  ColumnInitialize_double(column, 1000);
  fin = fopen(fname_out, "r");
  while(fscanf(fin, "%lf", &value)!=EOF){
    pos = (int)((value - shift)/binWidth);
    column[pos] = column[pos] + 1.0;
  }
  fclose(fin);

  for(i=0; i<1000; i++){
    column[i] = column[i]/dataNum;
  }

  fout = fopen(fname_hist, "w");
  for(i=0; i<1000; i++){
    if(column[i] > 0.0){
      fprintf(fout, "%lf\t%lf\n", shift + (i+0.5)*binWidth, column[i]);
    }
  }
  fclose(fout);

  plot_histogram_width_specified(fname_hist, "Elongation rate /min-1", "Frequency",
                 0.5, 1, 1, binWidth, 1.0);
  
}

void MeanVarAdd(double *mean, double *var, double value)
{
  (*mean)=(*mean)+value;
  (*var)=(*var)+pow(value, 2.0);
}

void MeanVarCalcFile(char fname[], double *mean, double *var)
{
  double value;
  double counter=0.0;
  FILE *fin=fopen(fname, "r");

  (*mean)=0.0;
  (*var)=0.0;

  while(fscanf(fin, "%lf", &value)!=EOF){
    MeanVarAdd(mean, var, value);
    counter=counter+1.0;
  }
  MeanVarFinal(mean, var, counter);
  fclose(fin);
}

void MeanVarFinal(double *mean, double *var, double counter)
{
  (*mean)=(*mean)/counter;
  (*var)=(*var)/counter-pow((*mean), 2.0);
}

void NormalizingFlValue(ROI roi[], int roiCount)
{
  char fname[] = "MeanFlForNormalization.dat";
  long lineNum;
  int i, j;
  int *sliceArray;
  double *meanFlTransition;
  FILE *fin;

  lineNum = FileLineNumberReturn(fname);
  meanFlTransition = (double *)malloc(sizeof(double)*lineNum);
  sliceArray = (int *)malloc(sizeof(int)*lineNum);
  if(meanFlTransition == NULL || sliceArray == NULL){
    printf("meanFlTransition pointer couldn't be assigned!\n");
    exit(1);
  }

  fin = fopen(fname, "r");
  for(i=0; i<lineNum; i++){
    fscanf(fin, "%d\t%lf", &sliceArray[i], &meanFlTransition[i]);
  }
  fclose(fin);

  for(i=0; i<roiCount; i++){
    if(roi[i].isConnected == 1){
      for(j=0; j<lineNum; j++){
    if(sliceArray[j] == roi[i].slice){
      roi[i].fl = roi[i].fl/meanFlTransition[j];
    }
      }
    }
  }

  free(meanFlTransition);
  free(sliceArray);
}

void NumberOrder(long number, int *order)
{
  if(number/10>=1){
    (*order)=(*order)+1;
    NumberOrder(number/10, order);
  }
}

void NumberStringReturn(char numberstring[], int rank, long number)
{
  int order=1;
  int i;
  char temp[100];

  NumberOrder(number, &order);//Determination of order
  if(rank<order){
    printf("Number is too large in 'NumberStringReturn'\n");
    exit(EXIT_FAILURE);
  }

  if(rank==order){
    sprintf(numberstring, "%ld", number);
  }else{
    for(i=0;i<rank-order;i++){
      temp[i]='0';
    }
    temp[rank-order]='\0';
    sprintf(numberstring, "%s%ld", temp, number);
  }
}

void PersistenceAnalysis(void)
{
  int choice;
  int loopindex = 0;

  while(loopindex == 0){
    printf("<PERSISTENCE ANALYSIS>\n");
    printf("1: Distribution of elongation rate for the last 10 time points\n");
    printf("2: Distribution of Fluorescence distribution at the last points\n");
    printf("3: Correlation between pre-exposure elongation rate and fluorescence intensity\n");
    printf("0: Exit\n");
    InputNonNegativeIntValue(&choice, "Enter number");
    switch(choice){
      case 0:
        loopindex = 1;
        break;
      case 1:
        LastTenPointElongationRateDist();
        break;
      case 2:
        LastFlDistribution();
        break;
      case 3:
        ElRateRpoSCorrelation();
        break;
      default:
        break;
    }
  }
}

void plot_histogram_width_specified(char filename[], char xlabel[], char ylabel[],
                    float fill_level, int eps_output_index, int datanum,
                    double binSize, double binRatio)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xtics %lf\n", binSize*5);
  fprintf(gp, "set mxtics 5\n");

  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);
  fprintf(gp, "set boxwidth %lf\n", binSize*binRatio); //Added

  for(i=0;i<datanum;i++){
    if(i==0){
      fprintf(gp, "plot '%s' using 1:2 w boxes fs solid %f notitle\n",
          filename, fill_level);
    }else{
      fprintf(gp, "replot '%s' using 1:%d w boxes fs solid %f notitle\n",
          filename, i+2, fill_level);
    }
  }

  //ChangeExtentionTo(fname_plt, ".plt");
  //fprintf(gp, "save '%s'\n", fname_plt);
  gpsave(gp, fname_plt);

  if(eps_output_index==1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void plot_linespoints(char filename[], char plottype[], char xlabel[], char ylabel[], int xlogscale, int ylogscale, int eps_output_index, int datanum)
{
  int i;
  char fname_plt[1000];
  FILE *gp;

  strcpy(fname_plt, filename);

  gp = popen(GNUPLOT_PATH, "w");
  if(gp == NULL){
    fprintf(stderr, "Oops, I can't find %s.", GNUPLOT_PATH);
    exit(EXIT_FAILURE);
  }

  gpStandardFormat(gp);
  fprintf(gp, "set xlabel '%s'\n", xlabel);
  fprintf(gp, "set ylabel '%s'\n", ylabel);

  if(xlogscale==1){
    fprintf(gp, "set logscale x\n");
  }
  if(ylogscale==1){
    fprintf(gp, "set logscale y\n");
  }

  if(strcmp(plottype, "l") == 0){
    for(i=0;i<datanum;i++){
      if(i==0){
    fprintf(gp, "plot '%s' using 1:2 w %s lw 2 notitle\n",
        filename, plottype);
      }else{
    fprintf(gp, "replot '%s' using 1:%d w %s lw 2 notitle\n",
        filename, i+2, plottype);
      }
    }
  }else if(strcmp(plottype, "p")==0){
    for(i=0;i<datanum;i++){
      if(i==0){
    fprintf(gp, "plot '%s' using 1:2 w %s ps 2 notitle\n",
        filename, plottype);
      }else{
    fprintf(gp, "replot '%s' using 1:%d w %s ps 2 notitle\n",
        filename, i+2, plottype);
      }
    }
  }else{
    for(i=0;i<datanum;i++){
      if(i==0){
    fprintf(gp, "plot '%s' using 1:2 w %s lw 2 ps 2 notitle\n",
        filename, plottype);
      }else{
    fprintf(gp, "replot '%s' using 1:%d w %s lw 2 ps 2 notitle\n",
        filename, i+2, plottype);
      }
    }
  }

  ChangeExtentionTo(fname_plt, ".plt");
  gpsave(gp, fname_plt);
  //fprintf(gp, "save '%s'\n", fname_plt);

  if(eps_output_index == 1){
    eps_output(gp, fname_plt);
  }

  fflush(gp);
  getchar();
  pclose(gp);
}

void RoiNextCellAssign(ROI roi[], int i)
{
  int j;

  j = roi[i].prevCell;

  if(j >= 0){
    if(roi[j].nextCellNum > 1){
      printf("More than 2 lineages are converging to roi[%d] in 'RoiNextCellAssign()'!\n", j);
      //printf("roi[%d].nextCellNum = %d\n", j, roi[j].nextCellNum);
      exit(EXIT_FAILURE);
    }else{
      roi[j].nextCell[roi[j].nextCellNum] = i;
      roi[j].nextCellNum = roi[j].nextCellNum + 1;
    }
  }
}

