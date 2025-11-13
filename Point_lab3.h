#include<iostream>
#include<vector>
#include<cmath>
#include<stdexcept>

using std::vector;

class Point {
private:
	vector<double> coordinates;

public:
	Point() {};
	Point(const vector<double>& coord) : coordinates(coord) {}

	Point(size_t n) : coordinates(n, 0.0) {};

	double dimension() const {
		return coordinates.size();
	}

	double& operator[](size_t index) {
		if (index >= dimension()) {
			throw std::out_of_range("Ind do not exist!");
		}
		return coordinates[index];
	}

	double norma() const {
		double sum = 0.0;
		for (double coord : coordinates) {
			sum += coord * coord;
		}
		return std::sqrt(sum);
	}

	Point operator+(const Point& other) {
		if (dimension() != other.dimension()) {
			throw std::invalid_argument("Ind do now exist!");
		}
		vector<double> result(coordinates.size());
		for (size_t i = 0; i < dimension(); i++) {
			result[i] = coordinates[i] + other.coordinates[i];
		}
		return Point(result);
	}

	Point operator-(const Point& other) const {
		if (dimension() != other.dimension()) {
			throw std::invalid_argument("Ind do not exist!");
		}
		vector<double> result(coordinates.size());
		for (size_t i = 0; i < dimension(); i++) {
			result[i] = coordinates[i] - other.coordinates[i];
		}
		return Point(result);
	}

	Point operator*(double scalar) const {
		vector<double> result(coordinates.size());;
		for (size_t i = 0; i < coordinates.size(); i++) {
			result[i] = coordinates[i] * scalar;
		}
		return Point(result);
	}

	friend Point operator*(double scalar, const Point& point) {
		return point * scalar;
	}
};
