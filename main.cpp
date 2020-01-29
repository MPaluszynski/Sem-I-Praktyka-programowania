#include <iostream>
#include <cstdlib>
#include <cmath>

#define MAX_N 10000000

int tab_int[MAX_N];
double tab_double[MAX_N];
int arrayOfPositions[MAX_N];
int q, valueOfSecondCounter;

using namespace std;

void read(int n)
{
    for(int i=0 ; i<n ; i++)
    {
        cin >> tab_int[i];
    }
    for(int i=0 ; i<n ; i++)
    {
        tab_double[i] = tab_int[i];
    }
}

void read2(double n)
{
    for(int i=0 ; i<n ; i++)
    {
        cin >> tab_double[i];
    }
}

int f0(int theArray[], int theArrayLength)
{
    int a, b, index = 0, counter = 0, secondCounter = 0;

    for(a = 0; a < theArrayLength ; a++)
    {
        if(theArray[a] <= theArray[index])
        {
            index = a;
            counter++;
        }
    };

    for(b = 0; b < theArrayLength ; b++)
    {
        if(theArray[b] == theArray[index])
        {
            arrayOfPositions[secondCounter] = b+1;
            secondCounter++;
        }
    };
    valueOfSecondCounter = secondCounter;

    return index;
};

void f1(int theArray[], int theArrayLength)
{
    int x, y, z, i, j, temp, minIndex, newLength = theArrayLength;

    int *newArray = new int[theArrayLength];
    int *sortedArray = new int[theArrayLength];

    for(x=0 ; x < theArrayLength ; x++)
    {
        newArray[x]=theArray[x];
    }

    for(y=0 ; y < theArrayLength ; y++)
    {
        minIndex = f0(newArray,newLength);
        sortedArray[y] = newArray[minIndex];
        for(z = minIndex ; z < newLength ; z++)
        {
            newArray[z]=newArray[z+1];
        }
        newLength--;
    }

    for(z=theArrayLength-1 ; z >= 0 ; z--)
    {
        cout<<sortedArray[z]<<" ";
    }
};

double f2(double theArray[], int theArrayLength, int flag)
{
    int i;
    double j=0.0, val;
    for(i=0 ; i<theArrayLength ; i++)
    {
        j += theArray[i]*theArray[i];
    }
    if(flag == 0)
    {
        val = floor(sqrt(j));
    }
    else if(flag == 1)
    {
        val = j;
    }
    return val;
};

double f3(double theArray[], int theArrayLength)
{
    int i;
    double sum=0.0, avg=0.0, A=0.0, temp=0.0, val;

    for(i=0 ; i<theArrayLength ; i++)
    {
        sum += theArray[i];
    }

    avg = (sum/theArrayLength);

    A = f2(theArray, theArrayLength, 1)/theArrayLength;

    val = floor(sqrt(A - (avg*avg)));

    return val;

};

void f4(int theArray[], int theArrayLength)
{
    int i, j, temp;
    j=theArrayLength-1;
    i=0;

    while(i<j)
    {
        temp=theArray[i];
        theArray[i]=theArray[j];
        theArray[j]=temp;
        i++;
        j--;
    }
    for(int i=0 ; i < theArrayLength ; i++)
    {
        cout<<theArray[i]<<" ";
    }
};

void f5(int theArray[], int theArrayLength)
{
    int i, j, n, isPrime = 1;

    for(i = 0 ; i<theArrayLength ; i++)
    {
        n = theArray[i];
        for(j = 2; j < n; ++j)
        {
            if((n%j) == 0)
            {
                isPrime=0;
                break;
            }
            else
            {
                isPrime=1;
            }
        }
        if (isPrime==1 && n != 1)
        {
            cout << "1" <<" ";
        }
        else
        {
            cout<<"0"<<" ";
        }
    }
}

void f6(int theArray[], int theArrayLength){
    long long int n1, n2, v=0, A;
    long long int *arrayOfX = new long long int[theArrayLength];
    long long int *arrayOfY = new long long int[theArrayLength];

    //6 10 91 84 -75 76 -88 67 -33 -49 54 17

    n1=theArrayLength/2;
    n2=theArrayLength;

    for(int x=0, i=0 ; x<n1, i < n2 ; x++, i += 2){

        arrayOfX[x] = theArray[i];
        //cout << "X" << x+1 << ": ";
        //cout << arrayOfX[x] << " ";
    }

    for(int j = 1, y=0 ; y<n1, j < n2 ; y++,j += 2){
        arrayOfY[y] = theArray[j];
        //cout << "Y" << y+1 << ": ";
        //cout << arrayOfY[y] << " ";
    }

    v=0;
    for(int k = 0 ; k < (n1-1) ; ++k){
        //cout << "arrayOfX " << k+1<< " = "  << arrayOfX[k] << endl;
        //cout << "arrayOfY " << k+2<< " = "  << arrayOfY[k+1] << endl;
        //cout << "arrayOfX " << k+2<< " = "  << arrayOfX[k+1] << endl;
        //cout << "arrayOfY " << k+1 << " = " << arrayOfY[k] << endl;
        //cout<<endl;
        //cout << arrayOfX[k] << " * " << arrayOfY[k+1] << " - " << arrayOfX[k+1] << " * " << arrayOfY[k]<<endl;
        v += ((arrayOfX[k]*arrayOfY[k+1])-(arrayOfX[k+1]*arrayOfY[k]));
    }
    v = v + ((arrayOfX[n1-1]*arrayOfY[0])-(arrayOfX[0]*arrayOfY[n1-1]));
    if(v<0){
            v = v *(-1);
        }
    A = (v)/2;
    cout << A << endl;
}

void f7(double theArray[], int theArrayLength)
{

    double a, b, c, d, delta, x1,x2,x3, f,g,h;
    a = theArray[0];
    b = theArray[1];
    c = theArray[2];
    d = theArray[3];

    if(a == 0)
    {
        //R-nie kwadratowe:
        delta = (c*c)-(4*b*d);

        if(delta>0)
        {
            x1 = (-c + sqrt(delta)) / (2*b);
            x2 = (-c - sqrt(delta)) / (2*b);

            if(x1 < x2)
                cout << floor(x1) << " " << floor(x2) << " " << endl;
            else
                cout << floor(x2) << " " << floor(x1) << " " << endl;
        }
        else if(delta == 0)
        {
            x1 = (-c + sqrt(delta)) / (2*b);
            cout<< floor(x1) << endl;
        }
        else
        {
            x1 = -c/(2*b);
            x2 = sqrt(-delta)/(2*b);
            cout << 0 << floor(x1) << "+" << floor(x2)<< "i" << " " << floor(x1)<< "-" << floor(x2) << "i" << " " << endl;
        }
    }
    else
    {
        //R-nie szescienne:
        f = ((c/a)-((b*b)/(3*a*a)));
        g = ((2*b*b*b)/(27*a*a*a))-(b*c)/(3*a*a)+(d/a);
        h = ((g*g)/4)+((f*f*f)/27);

        if(h>0)
        {
            x1 = cbrt((-g/2)+sqrt(h))+cbrt((-g/2)-sqrt(h))-(b/(3*a));
            cout<< x1 << endl;
        }
        else if (f==0 && g==0)
        {

            x1 = -cbrt(d/a);
            cout<< x1 << endl;
        }
        else
        {
            double i, j, k, m, n, p;
            i = sqrt(((g*g)/4)-h);
            j = cbrt(i);
            k = acos(-(g/(2*i)));
            m = cos(k/3);
            n = sqrt(3)*sin(k/3);
            p = (-b)/(3*a);

            x1 = 2*j*m+p;
            x2 = -j*(m+n)+p;
            x3 = -j*(m-n)+p;

            cout<< x1 << endl;
            cout<< x2 << endl;
            cout<< x3 << endl;
        }
    }
}

void f8(int theArray[])
{
    long long int n, v;
    n = theArray[0];
    v = (n*(n+1)*(n+2)*(3*n+5))/12;
    cout << v << endl;
}

void f9(int theArray[], int theArrayLength) {
    for (int i = 0; i < theArrayLength; i++) {

        unsigned int number = theArray[i];
        unsigned int count = 0;
        while (number) {
            count += number & 1;
            number >>= 1;
        }
        cout << count << " ";
    }
}

int main()
{
    int subprogram, n;

    while(cin >> subprogram >> n)
    {
        if(subprogram != 7)
        {
            read(n);
        }
        else if(subprogram == 7)
        {
            read2(n);
        }

        switch (subprogram)
        {
        case 0:
            f0(tab_int, n);
            for(q=0 ; q<valueOfSecondCounter ; q++)
            {
                cout<<arrayOfPositions[q]<<" ";
            }
            break;

        case 1:
            f1(tab_int, n);
            break;
        case 2:
            cout<< f2(tab_double, n, 0) << endl;
            break;
        case 3:
            cout<<f3(tab_double, n)<< endl;
            break;
        case 4:
            f4(tab_int, n);
            break;
        case 5:
            f5(tab_int, n);
            break;
        case 6:
            f6(tab_int, n);
            break;
        case 7:
            f7(tab_double, n);
            break;
        case 8:
            f8(tab_int);
            break;
        case 9:
            f9(tab_int, n);
            break;

        };
    };
    return 0;
}
