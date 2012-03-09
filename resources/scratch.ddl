create table files (id INTEGER PRIMARY KEY, fullname VARCHAR(1000), size INTEGER, last_used DATETIME, precious INTEGER, readers INTEGER, writers INTEGER, area INTEGER);
create table users (id INTEGER PRIMARY KEY, key VARCHAR(100), host VARCHAR(200), port INTEGER);
create table file_users (file_id INTEGER, user_id INTEGER, perm CHAR(2));
