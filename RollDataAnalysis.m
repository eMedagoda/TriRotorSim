clc
clear all
close all

data = importdata('../TriRotorDocuments/roll_test_data.csv');

figure(1)
plot(data(:,2),data(:,17))
hold on
plot(data(:,2),data(:,33),'r')

figure(2)
plot(data(:,2),data(:,22))
hold on
plot(data(:,2),data(:,23),'r')

figure(3)
plot(data(:,2),data(:,24))
hold on
plot(data(:,2),data(:,25),'r')

figure(4)
plot(data(:,2),data(:,20))
hold on

figure(5)
plot(data(:,2),data(:,27))
hold on
plot(data(:,2),data(:,28),'r')
