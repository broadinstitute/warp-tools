import pytest
from student import Student

@pytest.fixture
def student():
	s = Student("Alice", {"Math": 90, "English": 85, "History": 92})
	return s

def test_get_name(student):
	assert student.get_name() == "Alice"

def test_get_subjects(student):
	assert student.get_subjects() == ["Math", "English", "History"]

def test_get_grade(student):
	assert student.get_grade("Math") == 90
	assert student.get_grade("English") == 85
	assert student.get_grade("History") == 92

def test_set_name(student):
	student.set_name("Bob")
	assert student.get_name() == "Bob"

def test_set_grade(student):
	student.set_grade("Math", 95)
	assert student.get_grade("Math") == 95

def test_get_median_grade(student):
	assert student.get_median_grade() == 90

def test_add_subject(student):
	student.add_subject("Science", 87)
	assert student.get_subjects() == ["Math", "English", "History", "Science"]
	assert student.get_grade("Science") == 87

def test_get_lowest_grade_subject(student):
	assert student.get_lowest_grade_subject() == "English"

def test_get_highest_grade_subject(student):
	assert student.get_highest_grade_subject() == "History"

