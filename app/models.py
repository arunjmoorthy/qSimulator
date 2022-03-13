from operator import mod
from django.db import models

# Create your models here.


class Input(models.Model):
    random_paulis = models.IntegerField()
    physical_qubits = models.IntegerField()
    error_plot = models.BooleanField(null=True, blank=True, default=None)


class Graph(models.Model):
    number_of_cliques = models.IntegerField()
    text = models.TextField()
    image_one = models.ImageField(upload_to="graph")
    image_two = models.ImageField(upload_to="graph")
    text_one = models.TextField()
    text_two = models.TextField()

