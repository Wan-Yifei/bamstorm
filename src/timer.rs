use std::time::SystemTime;
use std::any::type_name;

pub fn type_name_of<T>(_: T) -> &'static str {
    type_name::<T>()
}

pub fn timeit<F: Fn() -> T, T>(f: F) -> T {
  let start = SystemTime::now();
  let result = f();
  let end = SystemTime::now();
  let duration = end.duration_since(start).unwrap();
  let func_name = type_name_of(f);
  println!(">>>> TIMER: Function {:?} executed in {:?} seconds:", func_name, duration.as_secs_f64());
  result
}
