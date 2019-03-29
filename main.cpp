#include<iostream>
#include<iomanip>

using namespace std;

class GaussianMatrix
{
    private:
        class GaussianMatrix *next_augmented;
        int rows;
        float *matrix;
        int state_level;
        bool solvable;
        float *solution;
    public:
        GaussianMatrix(float *mat, int r);
        void set_matrix(float *m);
        void print_matrix();
        bool augment();
        void solve();
        void print_solution();
};

GaussianMatrix::GaussianMatrix(float *mat, int r)
{
    matrix = (float *)malloc(sizeof(float)*r*(r+1));
    rows = r;
    next_augmented = NULL;
    state_level = -1;
    solvable = false;
    set_matrix(mat);
}

void GaussianMatrix::set_matrix(float *m)
{
    int i, elements;
    for(i = 0, elements = rows*(rows+1); i < elements; i++) { //use memcpy()
        *(matrix+i) = *(m+i);
    }
}

void GaussianMatrix::print_matrix()
{
    int i, elements;
    for(i = 0, elements = rows*(rows+1); i < elements; i++) {
        if(!(i%(rows+1))) {
            cout<<"\n";
        }
        cout<<left<<setw(12)<<*(matrix+i);
    }
    cout<<"\n##############################\n";
}

bool GaussianMatrix::augment()
{
    next_augmented = new GaussianMatrix(matrix, rows);
    class GaussianMatrix *temp = next_augmented;
    int pivot = state_level+1;
    int i, j, flag;
    float scalar, swap_temp;
    if(pivot == rows-1) {
        solvable = true;
        return true;
    }
    if(*(temp->matrix+pivot*(rows+1)+pivot)) {
        for(i = pivot+1; i < rows; i++) {
            scalar = *(temp->matrix+i*(rows+1)+pivot)/(*(temp->matrix+pivot*(rows+1)+pivot));
            for(j = pivot; j <= rows; j++) { //column is 1 greater than row
                (*(temp->matrix+i*(rows+1)+j)) -= scalar * (*(temp->matrix+pivot*(rows+1)+j));
            }
        }
        temp->state_level = state_level+1;
    } else {
        flag = 0;
        for(i = pivot+1; i < rows; i++) {
            if(*(temp->matrix+i*(rows+1)+pivot)) {
                for(j = pivot; j <= rows; j++) {
                    swap_temp = *(temp->matrix+i*(rows+1)+j);
                    *(temp->matrix+i*(rows+1)+j) = *(temp->matrix+pivot*(rows+1)+j);
                    *(temp->matrix+pivot*(rows+1)+j) = swap_temp;
                }
                flag = 1;
                break;
            }
        }
        if(!flag) {
            return false;
        }
        temp->state_level = state_level;
    }
    solvable = temp->augment();
    return solvable;
}

void GaussianMatrix::solve()
{
    GaussianMatrix *final_gauss_matrix;
    float *final_matrix;
    int i,j;
    float value;

    solvable = augment();

    if(solvable) {
        final_gauss_matrix = this;
        while(final_gauss_matrix->next_augmented) {
            final_gauss_matrix = final_gauss_matrix->next_augmented;
        }
        final_gauss_matrix->print_matrix();

        final_matrix = final_gauss_matrix->matrix;

        solution = (float *)malloc(rows*sizeof(float));
        for(i=rows-1; i>=0; i--) {
            value = 0;
            for(j=i+1; j < rows; j++) {
                value += (*(solution+j)) * (*(final_matrix+i*(rows+1)+j));
            }
            *(solution+i) = (*(final_matrix+i*(rows+1)+rows) - value) / (*(final_matrix+i*(rows+1)+i));
        }
    }
}

void GaussianMatrix::print_solution()
{
    int i;
    if(!solution) {
        cout<<"Cannot be solved"<<endl;
        return;
    }
    cout<<"Solution: (";
    for(i = 0; i < rows-1; i++) {
        cout<<*(solution+i)<<", ";
    }
    cout<<*(solution+i)<<")"<<"\n\n";
}

int main()
{
    int num, i , j;
    float *input_matrix;
    cout<<"Enter the number of equations:\n";
    cin>>num;
    input_matrix = (float *)malloc(sizeof(float)*num*(num+1));
    cout<<"Enter the coefficients [A1.x1+A2.x2+...+An.xn=C]:\n";
    for(i = 0; i < num; i++) {
        for(j = 0; j <= num; j++) { //column is 1 greater than row
            cin>>*(input_matrix+i*(num+1)+j);
        }
    }
    class GaussianMatrix gm(input_matrix, num);
    gm.print_matrix();
    gm.solve();
    cout<<"\n\n";
    gm.print_solution();
    cout<<"\n\n";
    return 0;
}