#include<iostream>
#include<vector>
#include<sstream>
#include<algorithm>
#include<complex>
#include<cmath>

/**
* @brief Класс полиномов.
*/
class Polynom
{
public:

	/**
	* @brief Конструктор класса.
	* @param number Число на вход.
	* @return Полином из числа.
	*/
	Polynom(unsigned int number);

	/**
	* @brief Метод преобразования в строку.
	* @return Полином в строку.
	*/
	std::string ToString();

	/**
	* @brief Геттер для степени полинма.
	* @return Степень полинома.
	*/
	size_t GetPower();

	/**
	* @brief Геттер для размера.
	* @return Размер полинома.
	*/
	size_t GetSize();

	/**
	* @brief Геттер для вектора коэффициентов.
	* @return Вектор коэффициентов.
	*/
	std::vector<unsigned int> GetKoef();

	/**
	* @brief Переопределения умножения.
	* @param other Полином.
	* @return Произведение полиномов.
	*/
	Polynom operator*(const Polynom& other);

	/**
	* @brief Переопределние сложения.
	* @param other Полином.
	* @return Сумма полиномов.
	*/
	Polynom operator+(const Polynom& other);

	/**
	* @brief Метод нахождения значения в точке.
	* @param dot Точка.
	* @return Значение в точке.
	*/
	unsigned int InDot(const unsigned int dot);

	/**
	* @brief Метод для переворота коэффициентов полинома.
	* @return Коэффициенты полинма в обратном порядке.
	*/
	std::vector<unsigned int> Reverse();

	/**
	* @brief Переопредление оператора вывода.
	* @param os Поток вывода.
	* @param polynom Полином.
	* @return Вывод в консоль.
	*/
	friend std::ostream& operator<<(std::ostream& os, Polynom& polynom);
private:
	unsigned int number;
	std::vector<unsigned int> koef;
};

Polynom::Polynom(unsigned int number) : number(number)
{
	size_t power = GetPower();

	for (size_t i = 0; i < power + 1; i++)
	{
		koef.push_back(number / static_cast<unsigned int>(std::pow(10, power - i)));
		number %= static_cast<unsigned int>(std::pow(10, power - i));
	}
}

std::string Polynom::ToString()
{
	std::stringstream buffer;

	for (size_t i = 0; i < koef.size(); i++)
	{
		buffer << koef[i] << "x^" << (koef.size() - i - 1);
		if (i != koef.size() - 1)
		{
			buffer << " + ";
		}
	}
	return buffer.str();
}

size_t Polynom::GetPower()
{
	std::ostringstream len;
	len << number;

	return len.str().size() - 1;
}

size_t Polynom::GetSize()
{
	return koef.size();
}

std::vector<unsigned int> Polynom::GetKoef()
{
	return koef;
}

Polynom Polynom::operator*(const Polynom& other)
{
	size_t resultPower = GetPower() * 2;
	Polynom result(0);

	result.koef.resize(resultPower + 1, 0);

	for (size_t i = 0; i < koef.size(); i++)
	{
		for (size_t j = 0; j < other.koef.size(); j++)
		{
			result.koef[i + j] += koef[i] * other.koef[j];
		}
	}

	return result;
}

Polynom Polynom::operator+(const Polynom& other)
{
	Polynom result(0);
	result.koef.resize(GetPower() + 1, 0);

	for (size_t i = 0; i < koef.size(); i++)
	{
		result.koef[i] = koef[i] + other.koef[i];
	}

	return result;
}

unsigned int Polynom::InDot(const unsigned int dot)
{
	unsigned int result = 0;

	for (size_t i = 0; i < koef.size(); i++)
	{
		result += koef[i] * pow(dot, koef.size() - i - 1);
	}

	return result;
}

std::vector<unsigned int> Polynom::Reverse()
{
	std::reverse(koef.begin(), koef.end());
	return koef;
}

std::ostream& operator<<(std::ostream& os, Polynom& polynom)
{
	os << polynom.ToString();
	return os;
}

int fact(int n)
{
	if (n < 0)
		return 0;
	if (n == 0)
		return 1;
	if (n > 0)
		return n * fact(n - 1);
}

size_t Delta(Polynom p, unsigned int k, size_t i)
{
	if ((k == 0))
		return p.InDot(0);
	if (i == 0)
		return p.InDot(k) - p.InDot(k - 1);

	return Delta(p, k, i - 1) - Delta(p, k - 1, i - 1);
}

void fft(std::vector<std::complex<double>>& poly, bool invert)
{
	int n = poly.size();
	if (n <= 1)
		return;
	std::vector<std::complex<double>> a0(n / 2), a1(n / 2);

	for (size_t i = 0; 2 * i < n; i++)
	{
		a0[i] = poly[2 * i];
		a1[i] = poly[2 * i + 1];
	}

	fft(a0, invert);
	fft(a1, invert);

	double ang = 2 * 3.14 / n * (invert ? -1 : 1);
	std::complex<double> w(1), wn(cos(ang), sin(ang));

	for (size_t i = 0; 2 * i < n; i++)
	{
		poly[i] = a0[i] + w * a1[i];
		poly[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
		{
			poly[i] /= 2;
			poly[i + n / 2] /= 2;
		}
		w *= wn;
	}
}

std::vector<unsigned int> IntToBinary(unsigned int num)
{
	std::vector<unsigned int> binary;

	if (num == 0)
	{
		binary.push_back(0);
	}
	else
	{
		while (num > 0)
		{
			binary.push_back(num % 2);
			num /= 2;
		}
	}

	reverse(binary.begin(), binary.end());

	return binary;
}

unsigned int BinaryToInt(const std::vector<unsigned int>& binary)
{
	int decimal = 0;
	int power = 1;

	for (int i = binary.size() - 1; i >= 0; i--)
	{
		if (binary[i] == 1)
		{
			decimal += power;
		}
		power *= 2;
	}

	return decimal;
}

std::vector<std::complex<double>> Multiply(const std::vector<std::complex<double>>& u, std::vector<std::complex<double>>& v)
{
	std::vector<std::complex<double>> fa(u.begin(), u.end()), fb(v.begin(), v.end());
	size_t n = 1;

	while (n < u.size() + v.size()) n <<= 1;
	fa.resize(n);
	fb.resize(n);

	fft(fa, false);
	fft(fb, false);

	for (size_t i = 0; i < n; i++)
	{
		fa[i] *= fb[i];
	}

	fft(fa, true);

	return fa;
}

std::vector<std::complex<double>> NumberToPoly(const std::vector<unsigned int>& num)
{
	std::vector<std::complex<double>> result(num.size());
	for (size_t i = 0; i < num.size(); ++i)
	{
		result[i] = num[num.size() - 1 - i];
	}
	return result;
}

std::vector<unsigned int> PolyToNumber(const std::vector<std::complex<double>>& poly)
{
	std::vector<unsigned int> result(poly.size());
	for (size_t i = 0; i < poly.size(); ++i)
	{
		result[i] = round(poly[i].real());
	}

	for (size_t i = 0; i < result.size(); ++i)
	{
		if (result[i] >= 2)
		{
			if (i + 1 == result.size()) result.push_back(0);
			result[i + 1] += result[i] / 2;
			result[i] %= 2;
		}
	}

	while (result.size() > 1 && result.back() == 0)
	{
		result.pop_back();
	}

	reverse(result.begin(), result.end());
	return result;
}

size_t methodKarazuba(const unsigned int a, const unsigned int b)
{
	Polynom pol1(a);
	Polynom pol2(b);

	Polynom result = pol1 * pol2;

	return result.InDot(10);
}

size_t methodToomKuk(const unsigned int u, const unsigned int v)
{
	Polynom ux(u);
	Polynom vx(v);
	Polynom wx = ux * vx;

	std::vector<unsigned int> res;
	res.resize(ux.GetSize() * 2, 0);

	for (size_t i = 0; i < res.size(); i++)
	{
		res[i] = Delta(wx, i, i - 1) / fact(i);
	}

	unsigned int result = 0;

	for (size_t i = 0; i < res.size(); i++)
	{
		unsigned int term = res[i];

		for (size_t j = 0; j < i; j++)
		{
			term *= (10 - j);
		}

		result += term;
	}

	return result;
}

size_t methodSchonhageStrassen(const unsigned int u, const unsigned int v)
{
	std::vector<unsigned int> ub = IntToBinary(u), vb = IntToBinary(v);
	std::vector<std::complex<double>> polyU = NumberToPoly(ub), polyV = NumberToPoly(vb);
	std::vector<std::complex<double>> resultPoly = Multiply(polyU, polyV);

	unsigned int result = BinaryToInt(PolyToNumber(resultPoly));

	return result;
}

int main()
{
	setlocale(LC_ALL, "Rus");

	int a = 26;
	int b = 89;

	std::cout << "Числа: " << a << " " << b << std::endl;
	std::cout << "Метод Карацубы: " << methodKarazuba(a, b) << std::endl;
	std::cout << "Метод Тоома-Кука: " << methodToomKuk(a, b) << std::endl;
	std::cout << "Метод Шенхаге-Штрассена: " << methodSchonhageStrassen(a, b) << std::endl;

	return 0;
}
