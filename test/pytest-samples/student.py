from __future__ import division

class Student:
	def __init__(self, name, grades_dict):
		self._name = name
		self._grades = grades_dict

	def get_name(self):
		return self._name

	def set_name(self, name):
		self._name = name

	def get_subjects(self):
		return list(self._grades.keys())

	def get_grade(self, subject):
		return self._grades.get(subject)

	def set_grade(self, subject, grade):
		if subject in self._grades:
			self._grades[subject] = grade

	def add_subject(self, subject, grade):
		self._grades[subject] = grade

	def get_lowest_grade_subject(self):
		return min(self._grades, key=self._grades.get)

	def get_highest_grade_subject(self):
		return max(self._grades, key=self._grades.get)
	
	# Function that computes the median grade of a student as a float
	def get_median_grade(self):
		grades = list(self._grades.values())
		grades.sort()
		if len(grades) % 2 == 0:
			return (grades[len(grades) // 2] + grades[len(grades) // 2 - 1]) / 2
		else:
			return grades[len(grades) // 2]
