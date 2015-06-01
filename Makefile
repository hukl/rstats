.PHONY: c compile

all:
	./rebar compile

c: compile
compile:
	./rebar compile skip_deps=true

shell: compile
	erl -pa ebin -pa deps/*/ebin

clean:
	rm -rf ebin/*.beam
