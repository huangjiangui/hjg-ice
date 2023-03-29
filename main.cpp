#include <iostream>
void simon(int n);
using namespace std;

int main(void)
{
	// using namespace std;
	simon(3);
	cout << "Pick an intege:"; 
	int n;
	cin >> n;
	simon(n);
	cout << "Done";
	return 0;
}
void simon(int n)
{

	cout << "simon says touch your toes" << n << " times." << endl;
}
