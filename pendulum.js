
let canvas = document.getElementById("canvas");
let ctx = canvas.getContext("2d");

let backgroundColor="rgb(0,0,0)";

var base_len;
var center;
var mass,length,theta,v,sum_mass;
var x,y,vx,vy;
var Energy;
var g=9.8;
var unit_step=500;
var count = 3;
var max_orbit_len = 40;
var show_orbit = true;

// var stroke_color;
var orbit_color;
var ball_radius;
var ball_color;
var rod_color;

var mode="standard"; //"standard" or other strings
var animID;

var start_T = null; 
var previous=null;
var start_E=null;

function randColor(){
  return `rgb(${parseInt(Math.random()*255)},${parseInt(Math.random()*255)},${parseInt(Math.random()*255)})`
}
var bugout = new debugout();


function init(){
  
  const base_angle = Math.PI/2;
  mass = Array(count);
  length = Array(count);
  theta = Array(count);
  v = Array(count);
  sum_mass = Array(count);
  ball_radius = Array(count);
  ball_color = Array(count);
  rod_color = Array(count);
  orbit_color = Array(count);
  x = Array(count);
  y = Array(count);
  vx= Array(count);
  vy= Array(count);

  past_orbit = NewArray(count);
  for(var i=0;i<count;i++){
    past_orbit[i] = new LoopQueue(max_orbit_len);
  }

  for(var i=0;i<count;i++){
    if(mode=="standard"){
    mass[i] = 1;
    length[i] = 1;
    theta[i]= base_angle*i;
    v[i] = 0;
    ball_radius[i] = 0.13*Math.sqrt(mass[i]);
    ball_color[i] = randColor();
    rod_color[i] = randColor();
    orbit_color[i] = ball_color[i];
    } else {
      mass[i] =1+ 2*Math.random();
      length[i] = 1+Math.random();
      theta[i] = 2*Math.PI *Math.random();
      v[i]=0;
      ball_radius[i] = 0.13*Math.sqrt(mass[i]);
      ball_color[i] = randColor();
      rod_color[i] = randColor();
      orbit_color[i] = ball_color[i];
    }
  }
  for(var i=count-1;i>=0;i--){
    if(i==count-1){
      sum_mass[i]=mass[i];
    }else {
      sum_mass[i]=sum_mass[i+1]+mass[i];
    }
  }
  get_xyvxvy();
}

function get_xyvxvy(){
  var tmpx=0;
  var tmpy=0;
  var tmpvx = 0;
  var tmpvy = 0;
  for(var i = 0; i < count;i++){
    tmpx += length[i]*Math.sin(theta[i]);
    tmpy += length[i]*Math.cos(theta[i]);
    tmpvx += length[i]*Math.cos(theta[i])*v[i];
    tmpvy += -length[i]*Math.sin(theta[i])*v[i];
    x[i]=tmpx;
    y[i]=tmpy;
    vx[i]=tmpvx;
    vy[i]=tmpvy;
  }
}

function getEnergy(){
  res= 0;
  for(var i=0;i<count;i++){
    res+= 1/2*mass[i]*(vx[i]**2+vy[i]**2);
    res -= mass[i]*g*y[i];
  }
  return res;
}
function getAcc(theta_,v_){ // 状态参量theta,v用参数传递，其他的用全局变量
  let M = NewMatrix(count);
  let b = NewArray(count);
  for(var i=0;i<count;i++){
    for(var j=0;j<count;j++){
      max_ = i>j?i:j;
      M[i][j]=length[i]*length[j]* Math.cos(theta_[i]-theta_[j])*sum_mass[max_];
    }
  }

  for(var i=count-1;i>=0;i--){
    b[i]=-sum_mass[i]*g*length[i]*Math.sin(theta_[i]);
  }
  for(var i=0;i<count;i++){
    for(var k=0;k<count;k++){
      max_ = i>k?i:k;
      b[i]-=sum_mass[max_]*length[i]*length[k]*v_[k]**2*Math.sin(theta_[i]-theta_[k]);
    }
  }
  return solve(M,b);// 解矩阵方程Mx = b
}
 // dtime以秒为单位计数
// 第一种方式
function evolve_1(dtime){
  acc = getAcc(theta,v);
  for(var i=0;i<count;i++){
    theta[i]+= v[i]*dtime;
    v[i]+= acc[i]*dtime ; 
  }
}
// 第二种方式
// function evolve_2(dtime){
//   acc = getAcc(theta,v);
//   for(var i=0;i<count;i++){
//     theta[i]+= (v[i]+acc[i]*dtime/2)*dtime;
//     v[i]+= acc[i]*dtime ; 
//   }
// }
// 第三种方式：四阶龙格库塔法
function evolve_rk4(dtime){
  // var acc1=NewArray(count),acc2=NewArray(count),acc3=NewArray(count),acc4=NewArray(count);
  var v1=NewArray(count),v2=NewArray(count),v3=NewArray(count);
  var theta1=NewArray(count),theta2=NewArray(count),theta3=NewArray(count);
  acc1 = getAcc(theta,v);
  for(var i=0;i<count;i++){
    v1[i]= v[i]+dtime/2*acc1[i];
    theta1[i] = theta[i]+dtime/2*v[i];
  }
  acc2 = getAcc(theta1,v1);
  for(var i=0;i<count; i++){
    v2[i] = v[i]+dtime/2*acc2[i];
    theta2[i] = theta[i]+dtime/2*v1[i];
  }
  acc3 = getAcc(theta2,v2);
  for(var i=0;i<count;i++){
    v3[i] = v[i]+dtime*acc3[i]
    theta3[i] = theta[i] + dtime* v2[i];
  }
  acc4 = getAcc(theta3,v3);
  for(var i=0;i<count;i++) {
    theta[i] += dtime/6*(v[i]+2*v1[i]+2*v2[i]+v3[i]);
    v[i] += dtime/6 * (acc1[i]+acc2[i]*2+acc3[i]*2+acc4[i]);
  }
}
unit_step=Math.round(unit_step/4); //减小计算量
evolve = evolve_rk4; //选用四阶rounge-kutta法
// bugout.log("method:rk_4")
// bugout.log(`unit_step:${unit_step}`)

function setSize() { //根据窗口大小重新调整
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
  min = Math.min(canvas.width,canvas.height)/2;
  base_len = min/(sum(length));
  center = Array(canvas.width/2,canvas.height/2);
}

function draw(){
  ctx.fillStyle = backgroundColor;
  ctx.fillRect(0,0,canvas.width,canvas.height);
  x_=center[0];
  y_=center[1];
  for(var i=0;i<count;i++){
    ctx.beginPath();
    if(i==0){
      ctx.moveTo(x_,y_);
    }  else{
      ctx.moveTo(x_+x[i-1]*base_len,y_+y[i-1]*base_len);
    }
    ctx.lineTo(x_+x[i]*base_len,y_+y[i]*base_len);
    ctx.lineWidth=4;
    ctx.strokeStyle=rod_color[i];
    ctx.stroke();
    ctx.closePath();
  }
  
  
  for(var i=0;i<count;i++){
    ctx.beginPath();
    ctx.arc(x_+x[i]*base_len,y_+y[i]*base_len,ball_radius[i]*base_len,0,2*Math.PI);
    ctx.fillStyle = ball_color[i];
    ctx.fill();
    ctx.closePath();
  }

  if(show_orbit){
    for(var i=0;i<count;i++){
      
      ctx.beginPath();
      flag=true;
      for(var j=0;j<max_orbit_len;j++){
        coordinate = past_orbit[i].get(j)
        if(coordinate!=undefined){
          // bugout.log(`${i},${j},${past_orbit[i].get(j)['x']}`)
          if(flag){
            flag=false;
            ctx.moveTo(x_+coordinate['x']*base_len,y_+coordinate['y']*base_len);
          } else {
            ctx.lineTo(x_+coordinate['x']*base_len,y_+coordinate['y']*base_len);
          }
        }
      }
      ctx.lineWidth=2;
      ctx.strokeStyle = orbit_color[i];
      ctx.stroke();
      ctx.closePath();
    }
  }

}


function update(timestamp) {
  // 以毫秒计数的时间戳
  if(previous==null){
    start_T=previous=timestamp;
    start_E=getEnergy();
    get_xyvxvy();

    draw();
  } else {
    if(timestamp-previous>100){
      // bugout.log(`WARNING:${timestamp-previous}`); 
      start_T += timestamp-previous; 
    } else{
      
      for(var i=0;i<unit_step;i++){
        evolve((timestamp-previous)/1000/unit_step);
      }
      get_xyvxvy();
      for(var i=0;i<count;i++){
        past_orbit[i].push({'x':x[i], 'y':y[i]});
      }
      Energy = getEnergy();
      
      // bugout.log(String((Energy-start_E)/start_E)+" "+String( timestamp-start_T));

    }
    previous=timestamp;
    draw();
  }
  animID = window.requestAnimationFrame(update);
}


window.onload = () => {
  setSize();
};

window.onresize = () => {
  setSize();
};
function remove(){
  delete mass,length,theta,v,sum_mass,ball_radius,ball_color,rod_color,orbit_color,x,y,vx,vy;
}
function start(){
  if(animID){
    remove();
    window.cancelAnimationFrame(animID);
  }
  init();
  setSize();
  window.requestAnimationFrame(update);
}
start();
