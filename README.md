# chemicalKinetic
Репозиторий для решения задачи химической кинетики.

chemkin_poly_first_order_pubo_and_qubo.py - решение PUBO и QUBO задачи для химической реакции первого порядка

chemkin_poly_second_order_pubo.py - решение PUBO задачи для химической реакции второго порядка

checmkin_poly_second_order_qubo.py - решение QUBO задачи для химической реакции второго порядка

chemkin_rms_second_order.py - МНК для химической реакции второго порядка

Задачи:

- [x] Лена:
  - [x] критерий сходимости Эйлера (2 и мб 2+)
  - [x] критерий обрезания ряда (2 и мб 2+)
- [x] Диана:
  - [x] презентация
  - [x] не забыть сделать акцент на особенностях данных
  - [x] не забыть про параметры anneal
- [x] Люба:
  - [x] графики при разных N_add для всех 4
  - [x] графики с шумом 1комп
  - [x] графики с 2 эксп точками
- [x] Маша:
  - [x] таблица по временам выполнения QUBO/PUBO
  - [x] найти ещё один параметр anneal
  - [x] cutoff для двух эксп точек
- [x] Денис:
  - [x] проверка кода
  - [x] рефакторинг кода
  - [x] подготовка репозитория
- [x] Миша:
  - [x] Обобщение QUBO на реакцию 2 компонент с распадами
  - [x] Обобщение QUBO/поиск механизмов обобщения для большего числа порядков
  - [x] График времени возведения k^m от m и от точности (длины вектора k)
- [ ] Сходить в бар  
  
- [x] запрогать, запустить и сравнить по времени toQubo
- [x] подумать, в каком виде вставлять графики результатов - мб в процентах
- [x] графики с шумом 2к

- [x] НЕ ЗАБЫТЬ что мы пубо к кубо сводим при сравнении времени выполнения
- [x] Потестить обрезание (в пубо? в кубо? в 1к? в 2к?)

- [x] Многомерный пубо (изи, но сделать нужно)
- [x] Не то что задача, подумать. Вот в двумерном кубо при увеличении N_add старшие степени к "вырождаются" - получаются колебания между двойкой и 0.125. При n_add = 7, 8 такое тоже происходило. Если будем N_add стремить к бесконечности - (возможно) такое тоже будет. Тогда мб имеет смысл среднее значение искать по старшим степеням к? Тогда получается концентрация, которая ближе к теоретической

