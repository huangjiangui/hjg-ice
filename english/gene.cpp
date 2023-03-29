#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>
using namespace std;

/*定义常量*/
#define SIZE 500 // 种群数量
#define MAX 1    // 个体最大值
#define MIN 2    // 个体最小值
/*定义运行参数*/
#define Cmax 100          // 确定函数最小值对应适应度函数的常数
#define Cmin 0            // 确定函数最大值对应适应度函数的常数
#define len1 9            // 决策变量x1字符串长度为9
#define len2 9            // 决策变量x1字符串长度为9
#define Chrom len1 + len2 //两个字符串整合成一个整体共18位，构成个体基因型
int FunctionMode = MAX;
int Size = 80;    // 初始化群体大小为80
int Max = 200;    //进化代数
double Pc = 0.6;  //设置交叉概率
double Pm = 0.01; //设置变异概率
/*个体数据结构*/
struct individual
{
  char chrom[Chrom + 1]; //一个个体染色体
  double value;          //个体的函数值
  double fitness;        //个体的适应度值
};
/*全局变量*/
int generation;                // 定义数据类型
int best_index;                // 最好个体下标
int worst_index;               // 最好个体下标
struct individual best;        // struct是一种结构体，可以包含多个不同类型的变量
struct individual worst;       //个体中最好的
struct individual currentbest; //当前最好的个体
struct individual pop[SIZE];
/*函数的原型*/
void initial();
void next();
void evaluate();
long chromosome(char *, int, int);
void value();
void fitness();
void BestandWrost();
void Perform();
void Select();
void Crossover();
void Mutation();
void Output();
int main()
{
  generation = 0; //初始化P(0)
  initial();      //产生初始群体，调用初始群体
  evaluate();     //初始群体适应度计算
  while (generation < Max)
  {
    generation++;
    next();
    evaluate();
    Perform();
    Output();
  }
  system("pause");
  return 0;
}
/*产生初代样本*/
void initial()
{
  int i, j;
  int flag = 1;
  srand((unsigned)time(NULL));
  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < Chrom; j++)
    {
      pop[i].chrom[j] = (rand() % 10 < 5) ? '0' : '1'; //随机数除以10取模运算
    }
    pop[i].chrom[Chrom] = '\0'; //字符串结束标志，防止乱码
  }
}
/*产生第一代*/
void next()
{
  Select();
  Crossover();
  Mutation();
}
/*评价个体*/
void evaluate()
{
  value();        // 计算目标函数值
  fitness();      // 计算适应度值
  BestandWrost(); // 找到最好和最差的个体
}
/*将二进制染色体解码为十进制整数.*/
long chromosome(char *string, int point, int length)
{
  int i;
  long decimal = 0L;
  char *pointer;
  for (i = 0, pointer = string + point; i < length; i++, pointer++)
  {
    decimal += (*pointer - '0') << (length - 1 - i);
  }
  return (decimal);
}
/*计算群体中每个个体的目标函数值*/
void value()
{
  int i;
  long temp1, temp2;
  double x1, x2;
  // Rosenbrock function
  for (i = 0; i < SIZE; i++)
  {
    temp1 = chromosome(pop[i].chrom, 0, len1); //二进制对应的十进制的值
    temp2 = chromosome(pop[i].chrom, len1, len2);
    x1 = 4.096 * temp1 / 1023.0 - 2.048;
    x2 = 4.096 * temp2 / 1023.0 - 2.048;
    pop[i].value = sin(x1 * x1) + (1 + x2) * (1 + x2) + exp(x1);
  }
}
/*前面已经计算好个体目标函数值，目标函数f(x)到适应度函数F(x)的转换关系*/
void fitness()
{
  int i;
  double temp;
  for (i = 0; i < SIZE; i++)
  {
    if (FunctionMode == MAX) //求目标函数的最大值，前面定义Cmin=0
    {
      if ((pop[i].value + Cmin) > 0.0)
      {
        temp = Cmin + pop[i].value;
      }
      else
      {
        temp = 0.0;
      }
    }
    else if (FunctionMode == MIN) //求目标函数最小值,Cmax=100
    {
      if (pop[i].value < Cmax)
      {
        temp = Cmax - pop[i].value;
      }
      else
      {
        temp = 0.0;
      }
    }
    pop[i].fitness = temp;
  }
}
void BestandWrost()
{
  int i;
  double sum = 0.0;
  // 找出当代最好和最差的个体
  best = pop[0];
  worst = pop[0];
  for (i = 1; i < SIZE; i++)
  {
    if (pop[i].fitness > best.fitness)
    {
      best = pop[i];
      best_index = i;
    }
    else if (pop[i].fitness < worst.fitness)
    {
      worst = pop[i];
      worst_index = i;
    }
    sum += pop[i].fitness;
  }
  // 找出最好的个体
  if (generation == 0)
  {
    currentbest = best;
  }
  else
  {
    if (best.fitness > currentbest.fitness)
    {
      currentbest = best;
    }
  }
}
/*用最好的取代最差的个体*/
void Perform()
{
  if (best.fitness > currentbest.fitness)
  {
    currentbest = pop[best_index];
  }
  else
  {
    pop[worst_index] = currentbest;
  }
}
/*基因选择*/
void Select()
{
  int i, index;
  double p, sum = 0.0;
  double cfitness[SIZE];
  struct individual newpopulation[SIZE];
  // 计算相对适应度
  for (i = 0; i < SIZE; i++)
  {
    sum += pop[i].fitness;
  }
  for (i = 0; i < SIZE; i++)
  {
    cfitness[i] = pop[i].fitness / sum;
  }
  //累计适应度
  for (i = 1; i < SIZE; i++)
  {
    cfitness[i] = cfitness[i - 1] + cfitness[i];
  }
  // 选择操作
  for (i = 0; i < SIZE; i++)
  {
    p = rand() % 1000 / 1000.0; //随机生成一个随机数除以1000得到0-999的余数再除以1000得到[0,0.999]的概率数
    index = 0;
    while (p > cfitness[index])
    {
      index++;
    }
    newpopulation[i] = pop[index];
  }
  for (i = 0; i < SIZE; i++)
  {
    pop[i] = newpopulation[i];
  }
}
/*基因交叉*/
void Crossover()
{
  int i, j;
  int index[SIZE];
  int point, temp;
  double p;
  char ch;
  //随机配对两个个体
  for (i = 0; i < SIZE; i++)
  {
    index[i] = i;
  }
  for (i = 0; i < SIZE; i++)
  {
    point = rand() % (SIZE - i);
    temp = index[i];
    index[i] = index[i + point];
    index[point + i] = temp;
  }
  //随机生成一个交叉点
  for (i = 0; i < SIZE - 1; i += 2)
  {
    p = rand() % 1000 / 1000.0;
    if (p < Pc)
    {
      point = rand() % (Chrom - 1) + 1;
      for (j = point; j < Chrom; j++)
      {
        ch = pop[index[i]].chrom[j];
        pop[index[i]].chrom[j] = pop[index[i + 1]].chrom[j];
        pop[index[i + 1]].chrom[j] = ch;
      }
    }
  }
}
/*基因变异*/
void Mutation()
{
  int i, j;
  double p;
  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < Chrom; j++)
    {
      p = rand() % 1000 / 1000.0;
      if (p < Pm)
      {
        pop[i].chrom[j] = (pop[i].chrom[j] == '0') ? '1' : '0';
      }
    }
  }
}
/*输出结果*/
void Output()
{
  int i;
  double sum;
  double average;
  //计算平均值
  sum = 0.0;
  for (i = 0; i < SIZE; i++)
  {
    sum += pop[i].value;
  }
  average = sum / SIZE;
  // print results of this pop
  printf("gen=%d,avg=%f,best=%f,", generation, average, currentbest.value);
  printf("chromosome=");
  for (i = 0; i < Chrom; i++)
  {
    printf("%c", currentbest.chrom[i]);
  }
  printf("\n");
}